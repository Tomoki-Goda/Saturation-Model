#ifndef R_FORMULA_HH
#define R_FORMULA_HH 1
#include<cmath>
#include<iostream>
#include<vector>
#include<string>
#include"control-default.h"
#include"constants.h"
#include"clenshaw.hh"
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
#if MODEL==1
template <typename ColG> class Sigma{
#else 
class Sigma{
#endif
// par and indivisual parameters are redundant. 
		//Collinear_Gluon xgpdf;
#if MODEL==1
		ColG xgpdf;
#endif
		double x2;
		double sigma_0,lambda, x_0, A_g,lambda_g,C,mu02,mu102,thresh_power;
		const double *par;
		inline double alpha(double mu2 ){
			static double b0= ((double)(33 -2*NF))/(12*PI);
			return( 1/(b0* log(mu2/LQCD2)));//LQCD2 lambda_QCD ^2
		}
		double Qs2(const double x,const double r);
	public:

		Sigma& operator=(const Sigma& rhs){
			init(rhs.par);
			return *this;
		} 
		explicit Sigma(void){ 
		}
		~Sigma(){
		}
		void set_x(const double &x);
		void init(const double * const &sigpar);
		double operator()(const double x, const double r) ;
};
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#endif

