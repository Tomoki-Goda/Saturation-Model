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
#if HANKEL==0
			const double r=rho/(1-rho);
#else
			const double r =rho;
#endif		
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
#if HANKEL==0
			return(val/pow(1-rho,2));
#else
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
double func_cheb(const double rho, void*p){
	const double r=rho/(1-rho);
	Sigma* sigma=(Sigma*)p;		
	double val=(*sigma)(r);///pow(1-rho,2);
	//printf("func %.3e r=%.3e, rho=%.3e\n",val,r, rho);getchar();
	return val;
}
double func_deriv(const double rho, void*p){
	gsl_function F;
 	F.function = func_cheb;
 	F.params=p;
	//const double r=rho/(1-rho);
	double val,err;
	gsl_deriv_forward(&F,rho,rho/10, &val, &err);
	
	return val;
}

class Laplacian_Sigma_Cheb{

	private:
		
		Sigma sigma;
		gsl_cheb_series *coeffs,*deriv;
		int r_npts=0;
		double x=0;
		char mode='l';
		void free_approx(){
			gsl_cheb_free (coeffs);
			gsl_cheb_free (deriv);
		}
		
	public:
		//double max=R_MAX, min=R_MIN;
		inline int set_kinem(double x){
			approximate(x);
			return 0;
		}	
		int approximate(const double x){
			this->x=x;
			sigma.set_kinem(x);
			//printf("Approximate\n");
			//getchar();
			gsl_function F;
 			F.function = func_cheb;
 			//F.function = func_deriv;
 			F.params=(void*)&sigma;
			//gsl_cheb_init(coeffs,&F,R_MIN,R_MAX);
			gsl_cheb_init(coeffs,&F,R_MIN/(1+R_MIN),R_MAX/(1+R_MAX));
			gsl_cheb_calc_deriv(deriv,coeffs);
			//printf("End\n");
			//getchar();
			return(0);
		}
		explicit Laplacian_Sigma_Cheb(){
			r_npts=0;
			
			//printf("sigma approx\n");
		}
		~Laplacian_Sigma_Cheb(){
			//printf("sigma approx end\n");
			free_approx();
			
		}
		void init(const int npts1,const double *par ,char mode){
			this->mode=mode;
			
			if(r_npts!=0){
				gsl_cheb_free(coeffs);
				gsl_cheb_free(deriv);
			}//else{
				//printf("first\n");
			r_npts=npts1;
			coeffs=gsl_cheb_alloc(npts1);
			deriv=gsl_cheb_alloc(npts1);
			sigma.init(par);
		}
		double operator()(const double rho)const{
			//const double r=rho/(1-rho);
			return(gsl_cheb_eval(deriv,rho));///pow(1-rho,2));
		}
		double operator()(const double rho, const std::vector<double> &par)const {
			const double r=rho/(1-rho);
			
			const double kt2=par[1];
			double val = 0;

//#if (LAPLACIAN==1||R_FORMULA==1)
			switch(mode){
				case 'l':
#if (IBP>=1 && HANKEL!=1)
#if IBP==1
					val=gsl_cheb_eval(deriv,rho);
					//val=gsl_cheb_eval(coeffs,rho);
					val*= sqrt(kt2)*r*std::cyl_bessel_j(1,r*sqrt(kt2));
#elif IBP==2
					val=gsl_cheb_eval(coeffs,rho);
					//val=gsl_cheb_eval(coeffs,rho);
					val*= kt2*r*std::cyl_bessel_j(0,r*sqrt(kt2));
#endif
#elif HANKEL==1
					val=gsl_cheb_eval(deriv,r);;
#else
					val=0;
#endif			
					break;
				case 's':
					val=gsl_cheb_eval(coeffs,r);
					val*=r*std::cyl_bessel_j(0,r*sqrt(kt2));
					break;
				default:
					printf("unknown option in laplacian sigma\n");
			}
		
			if(not(std::isfinite(val))){
				printf("2: val=%.3e for r= %.3e\n",val,r);
			}
			
			if(not(std::isfinite(val))){
				printf("3: val=%.3e for r= %.3e\nparameter: ",val,r);
				for(int i=0;i<par.size();i++){
					printf(" %.3e\t",par[i]);
				}printf("\n");
			}
			return(val);
		}
		
		double constant(double rho , const std::vector<double> &par)const {
			const double r=rho/(1-rho);
			const double kt2=par[1];
			double val;
			val=gsl_cheb_eval(deriv,rho)*pow(1-rho,2);
			val*=r*std::cyl_bessel_j(0,r*sqrt(kt2));
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
		CCIntegral cc=CCprepare(128,"dipole",50);
		
		inline double alpha(double mu2 )const{
			static double b0= ((double)(33 -2*NF))/(12*PI);
			return( 1/(b0* log(mu2/LQCD2)));//LQCD2 lambda_QCD ^2
		}
	public: 
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
			//this->x=x;
			//this->kt2=kt2;
			static double x_prev=0;
			if (x_prev!=x){
				x_prev=x;
				integrand.approximate(x);
			}
			const std::vector<double> par={x,kt2};
			double val=0;
#if GBW_APPROX==1
			if(x>0.7){
				//printf("approx\n");
				double qs2=(4*PI*PI*alpha(par[4])*par[1]*pow(x,-par[2])*pow(1-x,5.6))/(3*par[0]);
				val=2*par[0]*kt2/qs2*exp(-kt2/qs2);
				return(3.0/(8*PI*PI)*val);
			}
#endif
#if ADD_END>=0
			val=dclenshaw< const LSigma , const std::vector<double>& >(cc,integrand,par,R_MIN/(1+R_MIN),R_MAX/(1+R_MAX),INT_PREC/10,INT_PREC/100);
#endif
#if (IBP>=1&&ADD_END!=0)			
			val+=integrand.constant(R_MAX/(1+R_MAX),par);
			val-=integrand.constant(R_MIN/(1+R_MIN),par);
#endif

			return (3.0/(8*PI*PI)*val);
		}
		
		
};


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
				for(int i=0;i<kt2_npts;++i){
					kt2=((double)i)/(kt2_npts-1);
					kt2=kt2min*pow(4*kt2max/kt2min,kt2)/2;
					kt2_array[i] = kt2;
					aF_array[i+ j*kt2_npts] = aF(x,kt2,0);
				}
			}
				
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
			approximate(kt2max);
		}
		void init(const int npts1, const int npts2, const int npts3, const double *par ){
			x_npts=npts1;
			kt2_npts=npts2;
			
			x_array=(double*)malloc(x_npts*sizeof(double));
			kt2_array=(double*)malloc(kt2_npts*sizeof(double));
			aF_array=(double*)malloc(x_npts*kt2_npts*sizeof(double));
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

class Hankel_aF{
		Laplacian_Sigma sigma;
		//Interpolation
		double max_prev=0;
		//Dipole_Gluon aF;
		int kt2_npts,x_npts;
		gsl_interp_accel *x_accel_ptr, *kt2_accel_ptr;
		gsl_spline2d *  spline_ptr;
		double *r_array,*sigma_array,*kt2_array,*x_array,*aF_array;
		//double kt2min=1.0e-15,kt2max=-1;
		//Hankel
		const double rmax=R_MAX;
		//const int n=1000;
		gsl_dht* trans;
		int approximate(){
			std::vector<double> par(2,0);
#if IBP==0
			gsl_dht_init(trans,0,rmax);
#else
			gsl_dht_init(trans,1,rmax);
#endif
			for(int i=0;i<kt2_npts;i++){
				r_array[i]=gsl_dht_x_sample(trans,i);
				kt2_array[i]=gsl_dht_k_sample(trans, i);
				printf("%.3e\n",r_array[i]);
			}
			//printf("start\n");
			for(int j=0;j<x_npts;j++){
				x_array[j]=pow(10,-10+10*((double)j/(x_npts+1) ));
				sigma.set_kinem(x_array[j]);
				//printf("%d\n",j);
				for(int i=0;i<kt2_npts;i++){
					par[1]=kt2_array[i];
					//printf("sigma(x= %.3e, r= %.3e)\n",x_array[j],r_array[i]);
					sigma_array[i]=sigma(r_array[i],par);
					//printf("sigma(x= %.3e, r= %.3e)= %.3e\n",x_array[j],r_array[i],sigma_array[i]);
				}
				//for(int i=0;i<kt2_npts;i++){
				//	printf("%.3e\t%.3e\n",kt2_array[i],aF_array[j*kt2_npts+i] );
				//}
				gsl_dht_apply(trans,sigma_array,aF_array+j*kt2_npts );
			}
			
			gsl_spline2d_init (spline_ptr,kt2_array, x_array, aF_array, kt2_npts, x_npts);
			//printf("Hankel done\n");
			//getchar();
			return 0;
		}
	public:
		Hankel_aF(){
		}
		~Hankel_aF(){
			free_approx();
		}
		int export_grid(FILE* file){
			double val;
			for(int i=0;i<x_npts;i++){
			for(int j=0;j<kt2_npts;j++){
				val=gsl_spline2d_eval_extrap(spline_ptr, sqrt(kt2_array[j]), x_array[i], kt2_accel_ptr, x_accel_ptr);
				fprintf(file,"%.5e\t%.5e\t%.5e\n",x_array[i],kt2_array[j],val);
			}
			}	
			return 0;	
		}
		void free_approx(){
			gsl_spline2d_free (spline_ptr);
			gsl_interp_accel_free (x_accel_ptr);
			gsl_interp_accel_free (kt2_accel_ptr);
			free(kt2_array);
			free(x_array);
			free(sigma_array);
			free(r_array);
			free(aF_array);
			gsl_dht_free(trans);
		}
		void set_max(double kt2max){
			//this->kt2max=kt2max;
			approximate();
		}
		void init(const int npts1,const int npts2,const int npts3,const double *par ){
			
			x_npts=npts1;
			kt2_npts=npts2*100;
			
			trans=gsl_dht_alloc(kt2_npts);
			r_array=(double*)malloc(kt2_npts*sizeof(double));
			sigma_array=(double*)malloc(kt2_npts*sizeof(double));
			x_array=(double*)malloc(x_npts*sizeof(double));
			kt2_array=(double*)malloc(kt2_npts*sizeof(double));
			aF_array=(double*)malloc(x_npts*kt2_npts*sizeof(double));
			x_accel_ptr = gsl_interp_accel_alloc ();
			kt2_accel_ptr = gsl_interp_accel_alloc ();
			spline_ptr = gsl_spline2d_alloc(gsl_interp2d_bicubic,kt2_npts, x_npts); 
		//	aF.init(N_APPROX/2+50,par);
			sigma.init(npts3,par,'l');
		}
		double operator()(const double x,const double kt2,const double mu2)const{			
			double val = 0;
			val=gsl_spline2d_eval_extrap(spline_ptr, sqrt(kt2), x, kt2_accel_ptr, x_accel_ptr);
		//	printf("aF(%.3e, %.3e)= %.3e\n",x,kt2,val);
		//	getchar();
			return(3.0/(8*PI*PI)*val);
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
