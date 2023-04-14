///////////////////////////////////////////////////////
// Gluon Integrand.
// Thread safe
//
// USAGE 
// Gluon_Integrand integ
// integ.init(sigpar, mode, sig)// double*sigpar, mode is 'l' or 's' , Sigma sig
// integ(r, par) // double r is integ variable ch.var. is possible with R_CHANGE_VAR, doube* par={kt2,x}
// integ.constant(rmax,par)  // additional term for IBP  
//
// DSIGMA define at the end
//
// ///////////////////////////////////////////////////// 
#ifndef GLUON_INTEGRAND_HH
#define GLUON_INTEGRAND_HH 1

#include<cmath>
#include<iostream>
#include<vector>
#include<string>
#include<complex>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_spline.h>
#include"control.h"
#include"control-default.h"
#include"constants.h"
#include"clenshaw.hh"
#include"miscellaneous.hh"

#include"r-formula.hh"
#if MODEL==0
typedef Sigma_GBW SIGMA;
#else
typedef Sigma_BGK SIGMA;
#endif

class Gluon_Integrand{
		SIGMA *sigma;
		char mode='l';//l or s
		double ns_pow=500;
		const double *fixx;
		
		double sudakov(double r,double Q2,double k2){
			if(Q2<=k2){
				return 0;
			}
			double val=Q2*r*r/1.26095;
			val=0.2*3.0/(4*PI)*pow(log(val),2);
			return val;
		}
	public:
		Gluon_Integrand(SIGMA& sig){
			sigma=&sig;
		}
		void init(const double * const &par ,char mode);
		int set_x(const double &x);

		double operator()(const double rho, const std::vector<double> &par);
		//double constant(double r , const std::vector<double> &par);
};

#if GLUON_INTEGRAND_HH==2//no longer supported
template <typename Sig>  class Laplacian_Sigma{
	private:
		Sig *sigma;
		//double max=R_MAX, min=R_MIN; 
		double *r_array=NULL,*sigma_array=NULL;
		gsl_interp_accel *  r_accel_ptr;
		gsl_spline *  spline_ptr;
		char mode='l';//l or s
		int r_npts=0;
		int counter=0;
		const double *fixx;
		double sigma_0=0;
		int alloc_flag=0;
		double ns_pow=500;
		double rmax=R_MAX;
		
		void free_approx();
		void allocate(int npts1);
		int approximate(const double x);
		
	public:
		explicit Laplacian_Sigma(Sig & sig){
			sigma=&sig;
			r_npts=0;
			sigma_array=NULL;
			r_array=NULL;
			//printf("sigma approx\n");
		}
		~Laplacian_Sigma(){
			free_approx();
		}
		//double max=R_MAX, min=R_MIN;
		int set_x(const double& x);
		
		void init(const int npts1,const double * const &par ,char mode);
		double operator()(const double rho)const;
		int export_grid(FILE* file);
		double operator()(const double rho, const std::vector<double> &par);
		//double constant(double r , const std::vector<double> &par)const ;
};
#endif
#endif
