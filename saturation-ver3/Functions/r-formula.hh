#ifndef R_FORMULA_HH
#define R_FORMULA_HH 1
#include<cmath>
#include<iostream>
#include<vector>
#include<string>
#include"control.h"
#include"control-default.h"
#include"constants.h"
#include"gluons.hh"
#include"miscellaneous.hh"
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//   GBW / BGK dipoles
//   
//   Usage::
//   sigma.init(sigpar) //double* sigpar;
//   sigma(x,r) // double x,r 
//   
//   set_x( x ) can be used if gluon is approximated by chebyshev
// SIGMA 
////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if SIGMA_APPROX<0
typedef Chebyshev1D_Collinear_Gluon ColGlu;
#else 
typedef Collinear_Gluon ColGlu;
#endif




class Sigma{
	//Abstract base class
	protected:
		double x2;
		double xmax=1;
		double sigma_0,mu102,thresh_power;
		const double *par;
		//virtual
		 inline double freeze_x(const double);
		virtual double Qs2(const double x,const double r)const=0;
	public:
		virtual void set_x(const double &x){}
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
		double S(const double x, const double r) ;
};

///////////////////////////////////
//
///////////////////////////////////
class Sigma_GBW:public Sigma{
		double Qs2(const double x,const double r)const;
		double lambda, x_0;
	public:

		explicit Sigma_GBW(void){ 
		}
		~Sigma_GBW(){
		}
		void init(const double * const &sigpar);
};

///////////////////////////////////
//
///////////////////////////////////
class Sigma_BGK:public  Sigma{
		ColGlu xgpdf;
		double Qs2(const double x,const double r)const;
		double A_g,lambda_g,C,mu02;

//#if FREEZE_QS2>=1
//		const double xmax=0.8;
//#elif FREEZE_QS2==0
//	 	const double xmax=1;
//#endif		
		
		
	public:
		explicit Sigma_BGK(void){
#if SIGMA_APPROX<0
			xgpdf.allocate(N_CHEB);
#endif
#if FREEZE_QS2>=1
		xmax=0.98;
#elif FREEZE_QS2==0
	 	xmax=1;
#endif	
		}
		~Sigma_BGK(){
		}
		void set_x(const double &x);
		void init(const double * const &sigpar);
};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////



#endif

