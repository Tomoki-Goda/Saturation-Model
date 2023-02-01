#include<cmath>
#include<iostream>
#include<fstream>
#include<vector>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include"clenshaw.hh"
#include"gauss.hh"
//#include"./gluons.hh"
extern double INT_PREC;
extern int N_APPROX;
#ifndef LAPLACIAN
	#define LAPLACIAN 0
#endif

//FOR APPROXIMATION AND DERIVATIVES		
//class Laplacian_Sigma:public Sigma{
class Laplacian_Sigma{
	private:
		Sigma sigma;
		int r_npts;
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
	public:
		inline int set_kinem(double x){
			return(approximate(x));
		}	
		int approximate(const double x){
			this->x=x;
			//printf("approximate(%.3e)\n",x);
			//static int counter=0;
			double r;
			for (int j = 0; j < r_npts; j++){
				r=((double)j)/(r_npts-1);
				r=exp(std::log(R_MIN/2) + r * std::log(4*R_MAX/R_MIN));
				r_array[j]=r;
				sigma_array[j] = sigma(r,x);
				if(!isfinite(sigma_array[j])){
					printf("can not approximate sigma=%.3e\n",sigma_array[j] );
				}
			}
			gsl_spline_init (spline_ptr, r_array, sigma_array, r_npts);
			///printf("Approximated: %d %d\n" ,r_npts,++counter);
			//printf("\033[1A\033[2K\r");	
			return(0);
		}
		explicit Laplacian_Sigma(){
			//printf("sigma approx\n");
		}
		~Laplacian_Sigma(){
			//printf("sigma approx end\n");
			free_approx();
		}
		void init(const int npts1,const double (&par)[] ,char mode){
			this->mode=mode;
			if(sigma_array!=NULL){
				printf("refresh\n");
				free_approx();
			}
			r_npts=npts1;
			r_array=(double*)calloc(r_npts,sizeof(double));
			sigma_array=(double*)calloc(r_npts,sizeof(double));
			r_accel_ptr = gsl_interp_accel_alloc ();
			spline_ptr = gsl_spline_alloc(gsl_interp_cspline, r_npts); // cubic spline
			sigma.init(par);
		}
		inline double operator()(double r)const{
			return(gsl_spline_eval(spline_ptr, r,r_accel_ptr));
		}
		int export_grid(FILE* file){
			for(int i=0;i<r_npts;i++){
				fprintf(file,"%.5e\t%.5e\t%.5e\n",x,r_array[i],sigma_array[i]);
			}	
		return 0;	
		}
		double operator()(const double r, const std::vector<double> &par)const {
			const double x=par[0],kt2=par[1];
			double val = 0;

//#if (LAPLACIAN==1||R_FORMULA==1)
			switch(mode){
				case 'l':
					val=gsl_spline_eval_deriv2(spline_ptr, r,r_accel_ptr);
					val+=gsl_spline_eval_deriv(spline_ptr, r,r_accel_ptr)/r;
					break;
				case 's':
					val=gsl_spline_eval(spline_ptr, r,r_accel_ptr);
					break;
				default:
					printf("unknown option in laplacian sigma\n");
			}
		
		//	val=gsl_spline_eval(spline_ptr, r,r_accel_ptr);
//#elif LAPLACIAN==0
		//	val=gsl_spline_eval_deriv2(spline_ptr, r,r_accel_ptr);
		//	val+=gsl_spline_eval_deriv(spline_ptr, r,r_accel_ptr)/r;
//#endif//LAPLACIAN R_FORMULA
			if(not(std::isfinite(val))){
				printf("2: val=%.3e for r= %.3e\n",val,r);
			}
			
			val*=r*std::cyl_bessel_j(0,r*sqrt(kt2));
			if(not(std::isfinite(val))){
				printf("3: val=%.3e for r= %.3e\nparameter: ",val,r);
				for(int i=0;i<par.size();i++){
					printf(" %.3e\t",par[i]);
				}printf("\n");
			}
			return(val);
		}
};

class Dipole_Gluon{
	Laplacian_Sigma integrand;
	CCIntegral cc=CCprepare(128,"dipole",100);
	
	public: 
		//void init(Laplacian_Sigma* integrand ){
		inline void init(const int n,const double(&par)[] ){
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
			double val=3.0/(8*PI*PI)*dclenshaw< const Laplacian_Sigma , const std::vector<double>& >(cc,integrand,par,R_MIN,R_MAX,INT_PREC/10,INT_PREC/100);
		//	double val=3.0/(8*PI*PI)*dgauss< const Laplacian_Sigma , const std::vector<double>& >(integrand,par,R_MIN,R_MAX,INT_PREC/10,INT_PREC/100);
			return val;
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
				x=pow(10,-15+15*((double)j)/(x_npts-1));
				x_array[j] = x;
				for(int i=0;i<kt2_npts;++i){
					kt2=((double)i)/(kt2_npts-1);
					//r=pow(10,std::log10(R_MIN/2) + r * std::log10(4*R_MAX/R_MIN));
					//kt2=exp(std::log(kt2min/2)+kt2*std::log(4*(this->kt2max)/kt2min));
					//kt2=pow(10,std::log10(kt2min/2)+kt2*std::log10(4*(this->kt2max)/kt2min));
					kt2=kt2min*pow(4*kt2max/kt2min,kt2)/2;
					kt2_array[i] = kt2;
					aF_array[i+ j*kt2_npts] = aF(x,kt2,0);
					//printf("%.3e\n", kt2);
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
		void init(const int npts1,const int npts2,const double (&par)[] ){
			x_npts=npts1;
			kt2_npts=npts2;

			x_array=(double*)malloc(x_npts*sizeof(double));
			kt2_array=(double*)malloc(kt2_npts*sizeof(double));
			aF_array=(double*)malloc(x_npts*kt2_npts*sizeof(double));
			x_accel_ptr = gsl_interp_accel_alloc ();
			kt2_accel_ptr = gsl_interp_accel_alloc ();
			spline_ptr = gsl_spline2d_alloc(gsl_interp2d_bicubic,kt2_npts, x_npts); 
			aF.init( N_APPROX ,par);
		}
		double operator()(const double x,const double kt2,const double mu2)const{			
			double val = 0;
			//if(kt2>kt2max||kt2<kt2min){
			//	printf("%.3e<%.3e<%.3e\n",kt2min,kt2,kt2max);
			//}
			//if(x<1.0e-15||x>1){
			//	printf("%.3e<%.3e<%.3e\n",1.0e-15,x,1.0);
			//}
			val=gsl_spline2d_eval(spline_ptr,kt2, x,kt2_accel_ptr, x_accel_ptr);
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
