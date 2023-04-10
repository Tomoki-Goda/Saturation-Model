#ifndef R_FORMULA_HH
#define R_FORMULA_HH 1
#include<cmath>
#include<iostream>
#include<vector>
#include<string>
#include"control.h"
#include"control-default.h"
#include"constants.h"
#include"clenshaw.hh"
#include"gluons.hh"
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//   GBW / BGK dipoles
//   
//   Usage::
//   Sigma sigma // for BGK use Sigma<  (Cillinear_Gluon or Interpolate_Collinear_Gluon) >
// !! init for Collinear Gluon has to be fixed...
//   sigma.init(sigpar) //double* sigpar;
//   sigma(x,r) // double x,r 
//
// SIGMA 
////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if SIGMA_APPROX<0
typedef Chebyshev1D_Collinear_Gluon ColGlu;
#else 
typedef Collinear_Gluon ColGlu;
#endif




class Sigma{
// par and indivisual parameters are redundant. 
		//Collinear_Gluon xgpdf;
	protected:
		double x2;
		double sigma_0,mu102,thresh_power;
		const double *par;
		inline double alpha(double mu2 ){
			const double b0= ((double)(33 -2*NF))/(12*PI);
			return( 1/(b0* log(mu2/LQCD2)));//LQCD2 lambda_QCD ^2
		}
		virtual double Qs2(const double x,const double r)=0;
	public:
		virtual void set_x(const double x){}
		Sigma& operator=(const Sigma& rhs){
			init(rhs.par);
			return *this;
		} 
		explicit Sigma(void){ 
		}
		~Sigma(){
		}
		
		void init(const double * const &sigpar);
		double operator()(const double x, const double r) ;
};
class Sigma_GBW:public Sigma{
		double Qs2(const double x,const double r);
		double lambda, x_0;
	public:
		Sigma_GBW& operator=(const Sigma_GBW& rhs){
			init(rhs.par);
			return *this;
		} 
		explicit Sigma_GBW(void){ 
		}
		~Sigma_GBW(){
		}
		void init(const double * const &sigpar);
};
class Sigma_BGK:public  Sigma{
		ColGlu xgpdf;
		double Qs2(const double x,const double r);
		double A_g,lambda_g,C,mu02;
		
	public:
		Sigma_BGK& operator=(const Sigma_BGK& rhs){
			init(rhs.par);
			return *this;
		} 
		explicit Sigma_BGK(void){ 
		}
		~Sigma_BGK(){
		}
		void set_x(const double &x);
		void init(const double * const &sigpar);
};
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#endif

