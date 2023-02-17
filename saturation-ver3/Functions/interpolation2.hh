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
#include<pthread.h>
//#include<sys/wait.h>
//#include<unistd.h>
//#include"gauss.hh"
//#include"./gluons.hh"
extern double INT_PREC;
extern int N_APPROX;



//FOR APPROXIMATION AND DERIVATIVES		
//class Laplacian_Sigma:public Sigma{
class Laplacian_Sigma;
typedef struct {int j; Laplacian_Sigma *ptr;} sigmaopt;
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
		
		static void* compute(void*opt){
			sigmaopt* param=(sigmaopt*)opt;
			Laplacian_Sigma *sigmaptr=param->ptr;
			const int j=param->j;
			(sigmaptr->sigma_array)[j] = (sigmaptr->sigma)((sigmaptr->r_array)[j]);
			if(!std::isfinite((sigmaptr->sigma_array)[j])){
				printf("can not approximate sigma=%.3e\n",sigmaptr->sigma_array[j] );
			}
			
			return NULL;
			
		}
	public:
		//double max=R_MAX, min=R_MIN;
		inline int set_kinem(double x){
			approximate(x);
			return 0;
		}	
		int approximate_thread(const double x){
			sigma.set_kinem(x);
			pthread_t thread[r_npts];
			sigmaopt args[r_npts];
			int i1;//,i2,i3,i4;
			for (int j = 0; j < r_npts; j++){
				args[j].j=j;
				args[j].ptr=this;
				i1=pthread_create(thread+j,NULL,compute,(void*)(args+j) );
			}
			for(int i=0;i<r_npts;++i){
					pthread_join(thread[i],NULL);
			}
			//printf("ready\n ");
			gsl_spline_init (spline_ptr, r_array, sigma_array, r_npts);
				
			return(0);
		}
		int approximate(const double x){
			sigma.set_kinem(x);
			for (int j = 0; j < r_npts; j++){
				sigma_array[j] = sigma(r_array[j]);
				if(!isfinite(sigma_array[j])){
					printf("can not approximate sigma=%.3e\n",sigma_array[j] );
				}
			}
			gsl_spline_init (spline_ptr, r_array, sigma_array, r_npts);
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
			
			double r;
			for (int j = 0; j < r_npts; j++){
				r=((double)j)/(r_npts-1);
				r=min*pow(4.0*max/min,r)/2.0;
				r_array[j]=r;
			}
		}
		double operator()(const double rho)const{
			printf("rformula\n");
			const double r=rho/(1-rho);
			return(gsl_spline_eval(spline_ptr, r,r_accel_ptr)/pow(1-rho,2));
		}
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
			const double kt2=par[0];
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
			const double kt=sqrt(par[0]);
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
		CCIntegral cc=CCprepare(128,"dipole",4,4);
		double x;	

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
		void set_x(double x){
			this->x=x;
			integrand.approximate(x);
			//integrand.approximate_thread(x);
		}
		double operator()(const double x,const double kt2,const double mu2){
			if (this->x!=x){
				this->x=x;
				integrand.approximate(x);
			}
			return (*this)(kt2,mu2);

		}
		double operator()(const double kt2,const double mu2){
			Kahn accum=Kahn_init(3);
			double rmax=R_MAX,rmin=R_MIN;
			const std::vector<double> par={kt2};
			double val=0;
			Kahn_clear(accum);
#if GBW_APPROX==1
			if(x>0.7){
				double qs2=(4*PI*PI*alpha(par[4])*par[1]*pow(x,-par[2])*pow(1-x,5.6))/(3*par[0]);
				val=2*par[0]*kt2/qs2*exp(-kt2/qs2);
				return(3.0/(8*PI*PI)*val);
			}
#endif
#if ADD_END>=0
			rmax=50.0;
#if MODEL==1
			rmax/=pow(1-x,2.5);
#endif
			if(rmax>R_MAX||!std::isfinite(rmax)){
				//printf("rmax %.3e reduced to %.3e\n",rmax,R_MAX );
				rmax=R_MAX;
			}
			double scale=(2*PI)/sqrt(kt2),imin=rmin;
			int j=(int)(rmax/scale);
			if(j<3){
				j=3;
			}
			
#if IBP==0||IBP==2 //actually important. better cancellation for convergence.
			double imax=PI/(sqrt(kt2)*4);
#elif IBP==1
			double imax=3*PI/(sqrt(kt2)*4);
#endif
			
			int flag=0;
			//double val_prev=0;
			for(int i=0;i<j;++i){
				imax=imax+scale;
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
//////////////////////////////////////////
class Approx_aF;
typedef struct {int i, j; Approx_aF* ptr; } parallel_arg;

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
			free(kt2_array);
			free(x_array);
			free(aF_array);
		//	munmap((void*)aF_array,kt2_npts*x_npts*sizeof(double));
		//	munmap((void*)kt2_array,kt2_npts*sizeof(double));
		}
		int approximate(const double kt2max){
			this->kt2max=kt2max;
			clock_t time=clock();
		//	printf("APPROXIMATE\n1.0e-15 < kt2 < %.3e \n1.0e-10<x<1\n",ktmax);
		//	printf("kt2: %d, x: %d \n",kt2_npts,x_npts);
			//if(this->kt2max>kt2max ||  kt2max<(this->kt2max/10)){
//			int thread;
//#pragma omp parallel firstprivate( aF)
//		{
			double kt2,x;
//#pragma omp for 
			for (int j = 0; j < x_npts; ++j){
				//SAME x=pow(2,-30+30*((double)j)/(x_npts-1));
				x=pow(10,-8+8*((double)j)/(x_npts-1));
				//x=pow(((double)j+1.0)/(x_npts),2);
				x_array[j] = x;
				aF.set_x(x);	
				if(j!=0){
					printf("\033[1A\033[2K\r");
				}
				printf("[ ");
//#pragma omp parallel 
//			{
//#pragma omp for 
				for(int i=0;i<kt2_npts;++i){
					//thread=omp_get_thread_num();
					//printf("%d ",thread);
					kt2=((double)i)/(kt2_npts-1);
					kt2=kt2min*pow(4*kt2max/kt2min,kt2)/2;
					kt2_array[i] = kt2;
					aF_array[i+ j*kt2_npts] = aF(kt2,0);
					if((i/4)*4==i){
						printf("=");
					}
				}

	//		}
				printf("\033[2K\r");
				
				printf(" approxed x=%.2e\n", x);

			}
//		}
			printf("\033[1A\033[2K Grid done\n");
				
			gsl_spline2d_init (spline_ptr,kt2_array, x_array, aF_array, kt2_npts, x_npts);
			//}
			time-=clock();
			
			printf("%.2e sec to approx\n",-((double)time/CLOCKS_PER_SEC) );
			return(0);
		}
		static void* compute(void* par){
			//printf("func\n");
			parallel_arg* param=(parallel_arg*)par;
			int i=param->i,j=param->j;
			Approx_aF* to=param->ptr;
			double kt2=((double)i)/(to->kt2_npts-1);
			kt2=to->kt2min*pow(4*(to->kt2max)/(to->kt2min),kt2)/2;
			(to->kt2_array)[i] = kt2;
			(to->aF_array)[i+ j*(to->kt2_npts)] = (to->aF)(kt2,0);
			//printf("func end\n");
			//pthread_exit(NULL);
			return NULL;
		}
		int approximate_thread(const double kt2max){
			this->kt2max=kt2max;
			clock_t time=clock();
			
#pragma omp parallel firstprivate( aF)
		{
//
			double kt2,x;
			pthread_t thread[kt2_npts];
			parallel_arg args[kt2_npts];
			int i1;//,i2,i3,i4;
#pragma omp for 			
			for (int j = 0; j < x_npts; ++j){
				x=pow(10,-8+8*((double)j)/(x_npts-1));
				x_array[j] = x;
				aF.set_x(x);	
				if(j!=0){
					printf("\033[1A\033[2K\r");
				}
				printf("[ ");
				
				
				for(int i=0;i<kt2_npts;++i){
					args[i].i=i;
					args[i].j=j;
					args[i].ptr=this;
					i1=pthread_create(thread+i,NULL,compute,(void*)(args+i) );
				}					
				for(int i=0;i<kt2_npts;++i){
					pthread_join(thread[i],NULL);
				}
				printf("\033[2K\r");
				
				printf(" approxed x=%.2e\n", x);

			}
}
			printf("\033[1A\033[2K Grid done\n");
				
			gsl_spline2d_init (spline_ptr,kt2_array, x_array, aF_array, kt2_npts, x_npts);
			//}
			time-=clock();
			
			printf("%.2e sec to approx\n",-((double)time/CLOCKS_PER_SEC) );
			return(0);
		}

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
			//approximate(kt2max);
			approximate_thread(kt2max);
		}
		void init(const int npts1, const int npts2, const int npts3, const double *par ){
			x_npts=npts1;
			kt2_npts=npts2;
			
			x_array=(double*)malloc(x_npts*sizeof(double));
			kt2_array=(double*)malloc(kt2_npts*sizeof(double));
			//kt2_array=(double*)mmap(NULL,kt2_npts*sizeof(double),PROT_READ|PROT_WRITE,MAP_SHARED|MAP_ANONYMOUS,-1,0);
			aF_array=(double*)malloc(x_npts*kt2_npts*sizeof(double));
			//aF_array=(double*)mmap(NULL,x_npts*kt2_npts*sizeof(double),PROT_READ|PROT_WRITE,MAP_SHARED|MAP_ANONYMOUS,-1,0);
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
