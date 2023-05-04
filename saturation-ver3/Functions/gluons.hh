//  Claculating of xgpdf(x,Q^2) by DGLAP evolution of initial condition 
//  xgpdf(x,Q^2_0) 
// Originally written by S. Sapeta 
// Modified by T. Goda
// 
// COLGLU Defined at the end
 
#ifndef GLUONS_HH
#define GLUONS_HH
#include<cmath>
//#include"control.h"
//#include"control-default.h"
#include"constants.h"
#include"clenshaw.hh"
#include"polygamma.hh"
#include<complex>
#include<gsl/gsl_sf_gamma.h>
#include"chebyshev.hh"
//#define N_CHEB 25

class Collinear_Gluon{
	CCIntegral cc=CCprepare(256,"gluon",1,3);
	private:
		const double	beta = 6.6;
		double 			dgammafbeta;
		const double	n_0 = 0.5;       /// Maximal singluraity of integrand 
		std::complex<double> gammatilde(const std::complex<double>& n)const;
		
	public:
		explicit Collinear_Gluon(const Collinear_Gluon& init){
			dgammafbeta=init.dgammafbeta;
		}
		explicit Collinear_Gluon(){
			dgammafbeta=gsl_sf_gamma(beta)/PI;
		}
		~Collinear_Gluon(){
		}
		double operator()(const double y,const std::vector<double> &par)const;
		//MAIN FUNCTION
		double operator()(const double x, const double QQ,const double A_g,const double l_g)const;
};
///////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// 2D approximation with chebyshev
// not frequently used
///////////////////////////////////////////////////////////////////////
class Chebyshev_Collinear_Gluon{
	private:
		Collinear_Gluon xg;
	  	cheby cheb[2];
		double A_g=0,l_g=0;
		double xmin,xmax,q2min,q2max;
		const double *fixx=NULL;
		
	public:
		Chebyshev_Collinear_Gluon(unsigned n){
			const unsigned deg[2]={n,n};
			cheb[0]=PrepareChebyshev(deg,2);
			cheb[1]=PrepareChebyshev(deg[1],1);
		}
		~Chebyshev_Collinear_Gluon(){
			FreeChebyshev(cheb[0]);
			FreeChebyshev(cheb[1]);
		}
		double operator()(double *arg,const Collinear_Gluon& xg)const;
		
		int init(double xmin,double xmax,double q2min, double q2max, double A_g,double l_g );
		
		void set_x(const double &x);
		
		double operator()(const double x,const double Q2);
};
//////////////////////////////////////////////////////////////////////
// 1D approximation with chebyshev
///////////////////////////////////////////////////////////////////////
class Chebyshev1D_Collinear_Gluon{
	private:
		Collinear_Gluon xg;
	  	cheby cheb[1];
		double A_g=0,l_g=0;
		double q2min=0,q2max=0;
		const double *fixx=NULL;
		//unsigned deg[1]={0};
		unsigned flag=0;
	public:
		explicit Chebyshev1D_Collinear_Gluon(){
		}
		void allocate(unsigned n){
			if(flag==1){
				printf("1D chebyshev already allocated\n");
			}else{ 
				const unsigned deg[1]={n};
				cheb[0]=PrepareChebyshev(deg,1);
				flag=1;
			}
		}
		
		~Chebyshev1D_Collinear_Gluon(){
			FreeChebyshev(cheb[0]);
		}
		double operator()(double *arg,const Collinear_Gluon& xg)const;
		
		int init(double q2min, double q2max, double A_g,double l_g );
		
		void set_x(const double &x);
		
		double operator()(const double x,const double Q2, double A_g,double l_g )const;
};

//////////////////////////////////////////
//////////////////////////////////////////

#endif
