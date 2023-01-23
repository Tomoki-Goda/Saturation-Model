#include<cmath>
#include<iostream>
#include<fstream>
#include<vector>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include"clenshaw.h"

double par[]={25.0, 0.3,3.0e-4};
double Q0=1;
#define R_MIN 1.0e-6
#define R_MAX 50
double INT_PREC=1.0e-5;
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
};

class Laplacian_Sigma{
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
		
		const double* sigpar;
		//double x=0;
		
		int approximate(const double x){
			double r;
			for (int j = 0; j < r_npts; j++){
				r=pow(10,-10+13*((double)j)/(r_npts-1));
				r_array[j] = r;
				sigma_array[j] = sigma(r,x);
			}
			
			r_accel_ptr = gsl_interp_accel_alloc ();
			spline_ptr = gsl_spline_alloc(gsl_interp_cspline, r_npts); // cubic spline
			gsl_spline_init (spline_ptr, r_array, sigma_array, r_npts);
			return(0);
		}
	public:	
		void init(int npts1,const double (&par)[] ){
			r_npts=npts1;
			r_array=(double*)malloc(r_npts*sizeof(double));
			sigma_array=(double*)malloc(r_npts*sizeof(double));
			sigma.init(par);
		}
	
	private:
		double integrand(const double r, const std::vector<double> &par){
			double x=par[0],kt2=par[1];
			static double x_prev=0;
			if (x_prev!=x){
				x=x_prev;
				approximate(x);
			}	
			double val = gsl_spline_eval_deriv2(spline_ptr, r,r_accel_ptr);
			val+=gsl_spline_eval_deriv(spline_ptr, r,r_accel_ptr)/r;
			val*=r*std::cyl_bessel_j(0,r*sqrt(kt2));
			return(val);
		}
};
class aF{
	 Laplacian_Sigma integrand;
	 
	 
	public: 
		void init(int n,double (&par)[]){
			integrand.init(n,par);
		}
		double operator()(const double x,const double kt2){
			//this->x=x;
			//this->kt2=kt2;
			std::vector<double> par={kt2,x};
			double val=3.0/(4*PI*PI)*dclenshaw< Laplacian_Sigma , const std::vector<double>& >(integrand,par,R_MIN,R_MAX,INT_PREC);
			//double val=3.0/(4*PI*PI)*dclenshaw< double , const std::vector<double>& >(this->integrand,par,R_MIN,R_MAX,INT_PREC);
			return val;
		}
};

int main(){
	FILE* file=fopen("./test.txt","w");
	double Q2=10;
	double x=0.01;
	double val,kt2=0;
	aF gbw;
	gbw.init(250,par);
	for(int i =0;i<500;i++){
		kt2=pow(10,-3+5*((double)i)/(500));
		val=gbw(kt2,x);
		fprintf(file,"%.3e\t%.3e\n",kt2,val);
	}
	
	fclose(file);
	
		
	
}
