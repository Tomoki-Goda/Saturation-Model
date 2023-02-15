#include<cmath>
#include<iostream>
#include<fstream>
#include<vector>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include<gsl/gsl_dht.h>
#include<gsl/gsl_deriv.h>
#include<gsl/gsl_chebyshev.h>
#include"clenshaw.hh"
#include"gauss.hh"
#include<sys/wait.h>
#include<sys/mman.h>
#include<unistd.h>
//#include<pthread.h>
//#include"gauss.hh"
//#include"./gluons.hh"
extern double INT_PREC;
extern int N_APPROX;



//FOR APPROXIMATION AND DERIVATIVES		
//class Laplacian_Sigma:public Sigma{
class Laplacian_Sigma{
	private:
		double max=R_MAX, min=R_MIN; 
		Sigma sigma;
		int r_npts=0;
		gsl_interp_accel *  r_accel_ptr;
		gsl_spline *  spline_ptr;
		char mode='l';//l or s
		double *r_array=NULL,*sigma_array=NULL;
		
		
		void free_approx(){
			gsl_spline_free (spline_ptr);
			gsl_interp_accel_free (r_accel_ptr);
			free(r_array);
			free(sigma_array);
		}
		
		int counter=0;
		double x;
		double sigma_0=0;
	public:
		//double max=R_MAX, min=R_MIN;
		inline int set_kinem(double x){
			approximate(x);
			return 0;
		}	
		int approximate(const double x){
			this->x=x;
			//max=R_MAX/pow(1-x,4);
			//max=R_MIN/pow(1-x,2);
			//printf("approximate(%.3e)\n",x);
			//static int counter=0;
			double r;
			for (int j = 0; j < r_npts; j++){
				r=((double)j)/(r_npts-1);
				//r=exp(std::log(R_MIN/2) + r * std::log(4*R_MAX/R_MIN));
				r=min*pow(4.0*max/min,r)/2.0;
				r_array[j]=r;
				sigma_array[j] = sigma(r,x);
				if(!isfinite(sigma_array[j])){
					printf("can not approximate sigma=%.3e\n",sigma_array[j] );
				}
			}
			gsl_spline_init (spline_ptr, r_array, sigma_array, r_npts);
			//printf("Approximated: %d %d\n" ,r_npts,++counter);
			//printf("\033[1A\033[2K\r");	
			return(0);
		}
		explicit Laplacian_Sigma(){
			r_npts=0;
			sigma_array=NULL;
			r_array=NULL;
			//printf("sigma approx\n");
		}
		~Laplacian_Sigma(){
			//printf("sigma approx end\n");
			free_approx();
			
		}
		void init(const int npts1,const double *par ,char mode){
			sigma_0=par[0];
			this->mode=mode;
			if(r_npts!=0){
				free_approx();
			}//else{
				//printf("first\n");
			//}
			r_npts=npts1;
			r_array=(double*)calloc(r_npts,sizeof(double));
			sigma_array=(double*)calloc(r_npts,sizeof(double));
			r_accel_ptr = gsl_interp_accel_alloc ();
			spline_ptr = gsl_spline_alloc(gsl_interp_cspline, r_npts); // cubic spline
			sigma.init(par);
		}
		/*double operator()(const double rho)const{
			printf("rformula\n");
			const double r=rho/(1-rho);
			return(gsl_spline_eval(spline_ptr, r,r_accel_ptr)/pow(1-rho,2));
		}*/
		//int export_grid(FILE* file, FILE* file2){
		int export_grid(FILE* file){
			for(int i=0;i<r_npts;i++){
				fprintf(file,"%.5e\t%.5e\t%.5e\n",x,r_array[i],sigma_array[i]);
			}		
		return 0;	
		}
		double operator()(const double rho, const std::vector<double> &par)const {
#if (R_CHANGE_VAR==1)
			const double r=rho/(1-rho);
#elif (HANKEL==1||R_CHANGE_VAR==0)
			const double r =rho;
#endif		
			if(r>2*R_MAX){
				printf("too large, out of range %.4e - %.4e = %.4e\n",R_MAX,r,R_MAX-r);
			}
			if(r<R_MIN/2){
				printf("too small, out of range %.4e - %.4e = %.4e , rho=%.3e,%d\n",r,R_MIN,r-R_MIN,rho,R_CHANGE_VAR);
			}
			const double kt2=par[1];
			double val = 0;

//#if (LAPLACIAN==1||R_FORMULA==1)
			switch(mode){
				case 'l':
#if (IBP>=1 && HANKEL!=1)
#if IBP==1
					val=gsl_spline_eval_deriv(spline_ptr, r,r_accel_ptr);
					val*= sqrt(kt2)*r*std::cyl_bessel_j(1,r*sqrt(kt2));
#elif IBP==2
					val=gsl_spline_eval(spline_ptr, r,r_accel_ptr);
					val*=-kt2*r*std::cyl_bessel_j(0,r*sqrt(kt2));
#endif
#elif HANKEL==1
#if IBP==0
					val=gsl_spline_eval_deriv2(spline_ptr, r,r_accel_ptr);
					val+=gsl_spline_eval_deriv(spline_ptr, r,r_accel_ptr)/r;
#elif IBP==1
					val=sqrt(kt2)*gsl_spline_eval_deriv(spline_ptr, r,r_accel_ptr);
#endif
#else
					val=gsl_spline_eval_deriv2(spline_ptr, r,r_accel_ptr);
					val+=gsl_spline_eval_deriv(spline_ptr, r,r_accel_ptr)/r;
					val*=r*std::cyl_bessel_j(0,r*sqrt(kt2));
#endif			
					break;
				case 's':
					val=gsl_spline_eval(spline_ptr, r,r_accel_ptr);
					val*=r*std::cyl_bessel_j(0,r*sqrt(kt2));
					break;
				default:
					printf("unknown option in laplacian sigma\n");
			}
			//printf("2: val=%.3e for r= %.3e\n",val,r);
			if(not(std::isfinite(val))){
				printf("2: val=%.3e for r= %.3e\n",val,r);
			}
			
			if(not(std::isfinite(val))){
				printf("3: val=%.3e for r= %.3e\nparameter: ",val,r);
				for(int i=0;i<par.size();i++){
					printf(" %.3e\t",par[i]);
				}printf("\n");
			}
#if (R_CHANGE_VAR==1)
			return(val/pow(1-rho,2));
#elif (HANKEL==1||R_CHANGE_VAR==0)
			return val;
#endif
		}
		double constant(double rho , const std::vector<double> &par)const {
			const double r=rho/(1-rho);
			const double kt=sqrt(par[1]);
			double val;
#if IBP==1
			val=gsl_spline_eval_deriv(spline_ptr, r,r_accel_ptr);
			val*=r*std::cyl_bessel_j(0,r*kt);
#elif IBP==2
			val=kt*std::cyl_bessel_j(1,r*kt)*(gsl_spline_eval(spline_ptr,r,r_accel_ptr));
			val+=std::cyl_bessel_j(0,r*kt)*gsl_spline_eval_deriv(spline_ptr,r,r_accel_ptr);
			val*=r;
#endif

		//	if(fabs(val)>1.0e-3){
		//		printf("x= %.2e boundar= %.2e\n",x,val);
		//	}
			return(val);
		}
};


class Dipole_Gluon{
//#if IBP==0
		typedef Laplacian_Sigma  LSigma;
//#elif IBP==1
//	typedef Laplacian_Sigma_Cheb LSigma ;
//#endif
		const double *par;
		LSigma integrand;
		CCIntegral cc=CCprepare(128,"dipole",1,5);
		

		inline double alpha(double mu2 )const{
			static double b0= ((double)(33 -2*NF))/(12*PI);
			return( 1/(b0* log(mu2/LQCD2)));//LQCD2 lambda_QCD ^2
		}
	public: 
		Dipole_Gluon(){
		}
		~Dipole_Gluon(){
			
		}
		//void init(Laplacian_Sigma* integrand ){
		inline void init(const int n,const double *par ){
			this->par=par;
#if LAPLACIAN==0
			integrand.init(n,par,'l');	
#elif LAPLACIAN==1
			integrand.init(n,par,'s');
#endif	
		}
		double operator()(const double x,const double kt2,const double mu2){
			Kahn accum=Kahn_init(3);
		
			double rmax=R_MAX,rmin=R_MIN;
			//this->x=x;
			//this->kt2=kt2;
			static double x_prev=0;
			if (x_prev!=x){
				x_prev=x;
				integrand.approximate(x);
			}
			const std::vector<double> par={x,kt2};
			double val=0;
			Kahn_clear(accum);
#if GBW_APPROX==1
			if(x>0.7){
				//printf("approx\n");
				double qs2=(4*PI*PI*alpha(par[4])*par[1]*pow(x,-par[2])*pow(1-x,5.6))/(3*par[0]);
				val=2*par[0]*kt2/qs2*exp(-kt2/qs2);
				return(3.0/(8*PI*PI)*val);
			}
#endif
#if ADD_END>=0
			/*int j=(int)(2*Pi*(sqrt(kt2)+1));
			
			//int j=3;*/
			rmax=50.0/pow(1-x,2.5);
			if(rmax>R_MAX||!std::isfinite(rmax)){
				//printf("rmax %.3e reduced to %.3e\n",rmax,R_MAX );
				rmax=R_MAX;
			}
			//if(rmax<1){
			//	rmax=1;
			//}
			int j=(int)(sqrt(kt2)/(2*PI)*rmax);
			if(j<5){j=5;}
			
			//double scale=(rmax-rmin)/j,imin=rmin,imax=rmin;
			double scale=(2*PI)/sqrt(kt2),imin=rmin,imax=rmin;
			int flag=0;
			//double val_prev=0;
			for(int i=0;i<j;++i){
				imax=imin+scale;
				if(imax>rmax){imax=rmax;};
				//printf("min, max, %.2e %.2e\n",imin,imax);
#if R_CHANGE_VAR==1
				val=dclenshaw<const LSigma,const std::vector<double>&>(cc,integrand,par,imin/(1+imin),imax/(1+imax),INT_PREC/10,INT_PREC/100);
#elif R_CHANGE_VAR==0
				//printf("min, max, %.2e %.2e\n",imin,imax);
				val=dclenshaw<const LSigma,const std::vector<double>&>(cc,integrand,par,imin,imax,INT_PREC/10,INT_PREC/100);
#endif
				if(fabs(val)<1.0e-10){
					++flag;
					if(flag>3){
						break;//it is likely beyond this will be trivial
					}
				}else{
					flag=0;
				}
				//val_prev=val;
				
				accum+=val;
				imin=imax;
			}
			//printf("%d:\tx=%.3e kt2=%.3e, min=%.3e, max=%.3e,  rmin=%.3e, rmax=%.3e, scale=%.3e, val= %.2e\n",
				//j,x,kt2, imin,imax,rmin,rmax,scale,val);
			val=Kahn_total(accum);
#endif
#if (IBP>=1&&ADD_END!=0)			
			val+=integrand.constant(R_MAX/(1+R_MAX),par);
			val-=integrand.constant(R_MIN/(1+R_MIN),par);
#endif
			Kahn_free(accum);
			return (3.0/(8*PI*PI)*val);
		}
		
		
};

///////////FOR PTHREAD//////////////////
/*
typedef struct{
	int *args ;
	double* kt2_array,*x_array,*aF_array;
	double kt2max,kt2min;
	int kt2_npts;
	Dipole_Gluon* aFptr;
} aF_arg;
 
void* approximate_loop_recursive(void* arg){
	aF_arg* arg0=(aF_arg*)arg;
	//printf("branch \n");
	int i=arg0->args[0], child=arg0->args[1], j=arg0->args[2];
	if(child==pow(2,7)){
		getchar();
	}
	//int k,stat;
	if(child==1){
		
		//getchar();
		double kt2=((double)i)/((arg0->kt2_npts)-1);
		kt2=arg0->kt2min*pow(4*(arg0->kt2max)/(arg0->kt2min),kt2)/2;	
		(arg0->kt2_array)[i] = kt2;
		printf("%.2e %.2e %d\n", arg0->x_array[j],kt2, i+ j*(arg0->kt2_npts));
		//(arg0->aF_array)[i+ j*(arg0->kt2_npts)] = (*(arg0->aFptr))(arg0->x_array[j],kt2,0);
		printf("position=%d child=%d \n",i,child);
		
	}else{	
		printf("split %d \n",i);
		pthread_t thread1,thread2;
		aF_arg arg1=*arg0,arg2=*arg0;
		int args1[]={i-(child/4), child/2,j },args2[]={i+((child+2)/4), child/2,j};
		arg1.args=args1;
		arg2.args=args2;
		int i1,i2;
	
		i1=pthread_create(&thread1,NULL,&approximate_loop_recursive,(void*)&arg1);
		i2=pthread_create(&thread2,NULL,&approximate_loop_recursive,(void*)&arg2);
		pthread_join(thread1,NULL);
		pthread_join(thread2,NULL);
				
	}
}
*/
//////////////////////////////////////////

class Approx_aF{
	private:
		
		double max_prev=0;
		Dipole_Gluon aF;
		int kt2_npts,x_npts;
		gsl_interp_accel *x_accel_ptr, *kt2_accel_ptr;
		gsl_spline2d *  spline_ptr;
		double *kt2_array,*x_array,*aF_array;
		
		double kt2min=1.0e-15,kt2max=-1;
		
		void free_approx(){
			gsl_spline2d_free (spline_ptr);
			gsl_interp_accel_free (x_accel_ptr);
			gsl_interp_accel_free (kt2_accel_ptr);
			//free(kt2_array);
			free(x_array);
			//free(aF_array);
			munmap((void*)aF_array,kt2_npts*x_npts*sizeof(double));
			munmap((void*)kt2_array,kt2_npts*sizeof(double));
		}
		int approximate(const double kt2max){
			this->kt2max=kt2max;
			clock_t time=clock();
		//	printf("APPROXIMATE\n1.0e-15 < kt2 < %.3e \n1.0e-10<x<1\n",ktmax);
		//	printf("kt2: %d, x: %d \n",kt2_npts,x_npts);
			//if(this->kt2max>kt2max ||  kt2max<(this->kt2max/10)){
			double kt2,x;
			for (int j = 0; j < x_npts; ++j){
				//SAME x=pow(2,-30+30*((double)j)/(x_npts-1));
				x=pow(10,-8+8*((double)j)/(x_npts-1));
				//x=pow(((double)j+1.0)/(x_npts),2);
				x_array[j] = x;
				
				if(j!=0){
					printf("\033[1A\033[2K\r");
				}
				printf("[ ");
				for(int i=0;i<kt2_npts;++i){
					kt2=((double)i)/(kt2_npts-1);
					kt2=kt2min*pow(4*kt2max/kt2min,kt2)/2;
					kt2_array[i] = kt2;
					aF_array[i+ j*kt2_npts] = aF(x,kt2,0);
					if( ((10*i)/kt2_npts)*(kt2_npts)==10*i){
						printf("=");
					}
				}
				printf("\033[2K\r");
				
				printf(" approxed x=%.2e\n", x);
			}
			printf("\033[1A\033[2K Grid done\n");
				
			gsl_spline2d_init (spline_ptr,kt2_array, x_array, aF_array, kt2_npts, x_npts);
			//}
			time-=clock();
			
			printf("%.2e sec to approx\n",-((double)time/CLOCKS_PER_SEC) );
			return(0);
		}
		int approximate_loop(const double kt2max){
			this->kt2max=kt2max;
			clock_t time=clock();
			double x;
			for (int j = 0; j < x_npts; ++j){
				//printf("[ ");
				//printf("start splitting\n");
				x=pow(10,-8+8*((double)j)/(x_npts-1));
				x_array[j] = x;
				
				approximate_loop_recursive(kt2_npts/2, kt2_npts,j);
				printf("\033[2\r approxed x=%.2e", x);
			}
			printf("\033[1A\033[2K Grid done\n");
			gsl_spline2d_init (spline_ptr,kt2_array, x_array, aF_array, kt2_npts, x_npts);
			return 0;
		}	
		int approximate_loop_recursive(const int i, const int child, const int j){
			int k,stat;
			if(child==1){
				//printf("position=%d child=%d \n",i,child);
				double kt2=((double)i)/(kt2_npts-1);
				kt2=kt2min*pow(4*kt2max/kt2min,kt2)/2;	
				kt2_array[i] = kt2;
				aF_array[i+ j*kt2_npts] = aF(x_array[j],kt2,0);
			}else if(child==-1){
				exit(1);
			}else{
				k=fork();
				if(k==0){//child
					approximate_loop_recursive(i-(child/4), child/2,j);
					exit(0);
				}else{//parent
					approximate_loop_recursive(i+((child+2)/4), child/2,j);
					k=wait(&stat);
					//printf("%d,%d\n",k,stat);				
				}
			}
			return 0;
		}
		/*
		int approximate_loop(const double kt2max){
			aF_arg arg0;
			arg0.kt2min=kt2min;
			arg0.kt2max=kt2max;
			arg0.aF_array=aF_array;
			arg0.kt2_array=kt2_array;
			arg0.x_array=x_array;
			arg0.kt2_npts=kt2_npts;
			arg0.aFptr=&aF;
			double x;
			int args[3]={kt2_npts/2, kt2_npts,0};
			arg0.args=args;
			printf("af(1,10)=%.3e\n ",(*arg0.aFptr)(1,10,0));
			for (int j = 0; j < x_npts; ++j){
				x=pow(10,-8+8*((double)j)/(x_npts-1));
				x_array[j] = x;
				args[2]=j;
				printf("start branching\n");
				approximate_loop_recursive((void*)&arg0);
				printf("end branching\n");
				getchar();
			}
			
			gsl_spline2d_init (spline_ptr,kt2_array, x_array, aF_array, kt2_npts, x_npts);
			return 0;
		}*/	
		/*void* approximate_loop_recursive(void* arg){
			int*args0=(int*)args;
			int i=args0[0], child=args0[1], j=args0[2];
			
			int k,stat;
			if(child==1){
				double kt2=((double)i)/(kt2_npts-1);
				kt2=kt2min*pow(4*kt2max/kt2min,kt2)/2;	
				kt2_array[i] = kt2;
				aF_array[i+ j*kt2_npts] = aF(x_array[j],kt2,0);
			}else{
				pthread_t thread1,thread2;
				int args1[]={i-(child/4), child/2,j },args2[]={i+((child+2)/4), child/2,j};
				int i1,i2;
				
				i1=pthread_create(&thread1,NULL,&approximate_loop_recursive,(void*)args1);
				i2=pthread_create(&thread2,NULL,&approximate_loop_recursive,(void*)args2);
				pthread_join(thread1,NULL);
				pthread_join(thread2,NULL);
				
			}
		}
		*/
		
	public:
		int export_grid(FILE*file)const{
			for(int j=0;j< x_npts;j++){
				for(int i=0;i<kt2_npts;i++){
					fprintf(file ,"%.10e\t%.10e\t%.10e\n",x_array[j],kt2_array[i],aF_array[i+j*kt2_npts]);
				}
			}	
			return 0;
		}

		Approx_aF(){
		}
		~Approx_aF(){
			free_approx();
		}
		void set_max(double kt2max){
			this->kt2max=kt2max;
			approximate_loop(kt2max);
			//approximate(kt2max);
		}
		void init(const int npts1, const int npts2, const int npts3, const double *par ){
			x_npts=npts1;
			kt2_npts=npts2;
			
			x_array=(double*)malloc(x_npts*sizeof(double));
			//kt2_array=(double*)malloc(kt2_npts*sizeof(double));
			kt2_array=(double*)mmap(NULL,kt2_npts*sizeof(double),PROT_READ|PROT_WRITE,MAP_SHARED|MAP_ANONYMOUS,-1,0);
			//aF_array=(double*)malloc(x_npts*kt2_npts*sizeof(double));
			aF_array=(double*)mmap(NULL,x_npts*kt2_npts*sizeof(double),PROT_READ|PROT_WRITE,MAP_SHARED|MAP_ANONYMOUS,-1,0);
			x_accel_ptr = gsl_interp_accel_alloc ();
			kt2_accel_ptr = gsl_interp_accel_alloc ();
			spline_ptr = gsl_spline2d_alloc(gsl_interp2d_bicubic,kt2_npts, x_npts); 
			aF.init(npts3,par);
		}
		double operator()(const double x,const double kt2,const double mu2)const{			
			double val = 0;
			
			val=gsl_spline2d_eval(spline_ptr,kt2, x,kt2_accel_ptr, x_accel_ptr);
			//if(x<=0.5){
			//	val*=exp(-pow(1-x,-2));
			//}
			return(val);
		}
};

/*
double par[]={25.0, 0.3,3.0e-4};
double Q0=1;
#define R_MIN 1.0e-6
#define R_MAX 50
double INT_PREC=1.0e-5;

int main(){
	FILE* file=fopen("./test.txt","w");
	double Q2=10;
	double x=0.01;
	double val,kt2=0;
	Approx_aF gbw;
	gbw.init(100,100,par);
	gbw.set_max(1.0e+3);
	//std::vector<double> param(2,0);
	//param[0]=x;
	//param[1]=1;

	for(int i =0;i<500;i++){
		kt2=pow(10,-4+7*((double)i)/(500));
		//val=gbw(kt2,param);
		val=gbw(x,kt2);
		fprintf(file,"%.3e\t%.3e\n",kt2,val);
		//printf("%.3e\t%.3e\n",kt2,val);
	}
	
	fclose(file);
	
		
	
}
*/
