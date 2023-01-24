#include<cmath>
#include<iostream>
#include<fstream>
#include<vector>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include"clenshaw.hh"
#include"./gluons.hh"
extern double INT_PREC;


class Sigma{
		const double* sigpar;
		double x=0;
	public:
		void init(const double (&par)[]){
			sigpar=par;
		}
		explicit Sigma(void){
		}
		~Sigma(){
		}
		double operator()(const double r,const double x){//,const double Q2,const double*sigpar)const {
		 	this->x=x;
		 	return((*this)(r));
		 }
#if MODEL==0
		double operator()(const double r)const{
			double sigma_0=sigpar[0];
			double lambda=sigpar[1];
			double x_0=sigpar[2];
			if(x_0<1.0e-5||x_0>1.0e-3){
				return 0;
			}
			if(lambda<0.05||lambda>0.95){
				return 0;
			}
			
			return( sigma_0*(1-exp( - pow(r * Q0, 2) * pow(x_0/x, lambda)/4)) );	
		}
#elif MODEL==1
	private:
		Collinear_Gluon xg;
		
		inline double alpha(double mu2 ){
			static double b0= ((double)(33 -2*NF))/(12*PI);
			return( 1/(b0* log(mu2/LQCD2)));//LQCD2 lambda_QCD ^2
	
		}
		//inline PREC alpha(const PREC mu2)const{
		//	return 4.0/(9.0 *log( ((mu2>2*LQCD2)?(mu2):(2.0*LQCD2))/LQCD2));
		//}
		double operator()(const double r){
			double sigma_0=sigpar[0];
			
			//double A_g=sigpar[1];
			//double lambda_g=sigpar[2];
			if(sigpar[1]<0){
				return 0;
			}
			xg.set_xg_parameter(sigpar[1],sigpar[2]);
			double C=sigpar[3];
			double mu02=sigpar[4];
		
			double mu2=C/(r*r)+mu02;
			if(mu2<1){
				return(0);
			}
			return( sigma_0*(1-exp(-pow(r*PI,2)*alpha(mu2)*xg(x,mu2 )/(3*sigma_0) ) )); 	
		}
#endif		
};
//FOR APPROXIMATION AND DERIVATIVES		
//class Laplacian_Sigma:public Sigma{
class Laplacian_Sigma{
	private:
		Sigma sigma;
		int r_npts;
		gsl_interp_accel *  r_accel_ptr;
		gsl_spline *  spline_ptr;
		double *r_array,*sigma_array;
		
		
		void free_approx(){
			gsl_spline_free (spline_ptr);
			gsl_interp_accel_free (r_accel_ptr);
			free(r_array);
			free(sigma_array);
		}
		
		int approximate(const double x){
			double r;
			for (int j = 0; j < r_npts; j++){
				r=((double)j)/(r_npts);
				r=exp(std::log(R_MIN/2) + r * std::log(4*R_MAX/R_MIN));
				r_array[j]=r;
				sigma_array[j] = sigma(r,x);
			}
			gsl_spline_init (spline_ptr, r_array, sigma_array, r_npts);
			return(0);
		}
	public:	
		Laplacian_Sigma(){
		}
		~Laplacian_Sigma(){
			free_approx();
		}
		void init(const int npts1,const double (&par)[] ){
			r_npts=npts1;
			r_array=(double*)malloc(r_npts*sizeof(double));
			sigma_array=(double*)malloc(r_npts*sizeof(double));
			r_accel_ptr = gsl_interp_accel_alloc ();
			spline_ptr = gsl_spline_alloc(gsl_interp_cspline, r_npts); // cubic spline
			sigma.init(par);
		}
		double operator()(const double r, const std::vector<double> &par){
			static double x_prev=0;
			const double x=par[0],kt2=par[1];
			//std::cout<<std::scientific<<x<<"  "<<kt2<<std::endl;
			//static double x_prev=0;
			if (x_prev!=x){
				//printf("approx x_prev= %.3e-> x= %.3e\n",x_prev,x);
				x_prev=x;
				approximate(x);
			}	
			double val = 0;
			//val+=gsl_spline_eval(spline_ptr, r,r_accel_ptr);
			val+=gsl_spline_eval_deriv2(spline_ptr, r,r_accel_ptr);
			val+=gsl_spline_eval_deriv(spline_ptr, r,r_accel_ptr)/r;
			val*=r*std::cyl_bessel_j(0,r*sqrt(kt2));
			//printf("%.3e\n", val);
			return(val);
		}
};

class Dipole_Gluon{
	Laplacian_Sigma integrand;
	public: 
		//void init(Laplacian_Sigma* integrand ){
		void init(const int n,const double(&par)[] ){
			integrand.init(n,par);	
		}
		double operator()(const double x,const double kt2,const double mu2){
			//this->x=x;
			//this->kt2=kt2;
			const std::vector<double> par={x,kt2};
			double val=3.0/(8*PI*PI)*dclenshaw< Laplacian_Sigma , const std::vector<double>& >(integrand,par,R_MIN,R_MAX,INT_PREC/10,INT_PREC/100);
			//double val=integrand(kt2,par);
			//double val=3.0/(4*PI*PI)*dclenshaw< double , const std::vector<double>& >(this->integrand,par,R_MIN,R_MAX,INT_PREC);
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
		clock_t time=clock();
		//	printf("APPROXIMATE\n1.0e-15 < kt2 < %.3e \n1.0e-10<x<1\n",ktmax);
		//	printf("kt2: %d, x: %d \n",kt2_npts,x_npts);
			//if(this->kt2max>kt2max ||  kt2max<(this->kt2max/10)){
			double kt2,x;
			for (int j = 0; j < x_npts; ++j){
				x=pow(10,-10+10*((double)j)/(x_npts-1));
				x_array[j] = x;
				for(int i=0;i<kt2_npts;++i){
					kt2=((double)i)/(kt2_npts);
					//r=pow(10,std::log10(R_MIN/2) + r * std::log10(4*R_MAX/R_MIN));
					kt2=exp(std::log(kt2min/2)+kt2*std::log(4*(this->kt2max)/kt2min));
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
			aF.init(200,par);
		}
		double operator()(const double x,const double kt2,const double mu2)const{			
			double val = 0;
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
