#include<iostream>
#include<fstream>
#include<cmath>

#include<fstream>
#include<vector>
#include<chrono>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include<gsl/gsl_dht.h>
#include<gsl/gsl_deriv.h>
#include<gsl/gsl_chebyshev.h>
//#include"Functions/r-formula.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sum.h>
#define N_CHEB 20
#include"Functions/gluons.hh"
#include"Functions/Kahn.hh"
#include"Functions/clenshaw.hh"
#include"Functions/control-default.h"
#include"Functions/constants.h"

class testfunc0{
	
	public:
		double operator()(double r, double* p){
			double qs2=pow(2.0e-4/p[1],0.2); 
			
			return(std::cyl_bessel_j(0,r *p[0] )* (1-exp(-r*r *qs2/4))/r);
		}

};

class testfunc{
	Collinear_Gluon xg;
	inline double alpha(double mu2 ){
	
			static double b0= ((double)(33 -2*NF))/(12*PI);
			return( 1/(b0* log(mu2/LQCD2)));//LQCD2 lambda_QCD ^2
		}
	
	public:
		double operator()(double r, double* p){
			double q2=pow(r,-2)+1;
			//printf("xg(%.3e %.3e)=%.3e\n", p[1],q2,xg(p[1],q2,1,0.2));
			double qs2=4*PI*PI*alpha(q2)*xg(p[1],q2,1.0,0.2)/(3*85.0); 
			
			return(std::cyl_bessel_j(0,r *p[0] )* (1-exp(-r*r *qs2/4))/r);
		}

};


class testfunc_cheb{
	Chebyshev1D_Collinear_Gluon *xg;
	
	inline double alpha(double mu2 ){
	
			static double b0= ((double)(33 -2*NF))/(12*PI);
			return( 1/(b0* log(mu2/LQCD2)));//LQCD2 lambda_QCD ^2
		}
	
	public:
		testfunc_cheb(Chebyshev1D_Collinear_Gluon& g){
			xg=&g;
		}
		double operator()(double r, double* p){
			double q2=pow(r,-2)+1;
			//printf("xg(%.3e %.3e)=%.3e\n", p[1],q2,xg(p[1],q2,1,0.2));
			double qs2=4*PI*PI*alpha(q2)*(*xg)(p[1],q2)/(3*85.0); 
			//printf("qs2=%.3e r=%.3e x=%.3e, k=%.3e\n", qs2,r,p[1],p[0]);
			return(std::cyl_bessel_j(0,r *p[0] )* (1-exp(-r*r *qs2/4))/r);
		}

};

int main(int argc, char** argv ){
	double par[10]={0};
	CCIntegral cc=CCprepare(64,"dipole",1,5);
	
	par[0]=atof(argv[1]);
	par[1]=atof(argv[2]);
	double imin,imax;
	Kahn acc=Kahn_init(3);
	double val=0,sum=0;
	double scale=1*(2*PI)/sqrt(par[0]);
	
	double total=0;
	Chebyshev1D_Collinear_Gluon xg;
	
	testfunc_cheb func1(xg);
	testfunc func2;
	testfunc0 func3;
	
	
	Kahn_clear(acc);
	sum=0;
	imin=1.0e-8;
	imax=10*PI/(4*par[0]);
	
	printf("\n***************************\nGBW\n****************************\n");
	for (int i=0;i<10;++i){
		imax+=scale;
		val=dclenshaw<testfunc0,double*>(cc, func3, par,imin,imax,1.0e-10,1.0e-10);
		sum+=val;
		acc+=val;
		total=Kahn_total(acc);
		printf("[%.3e %.3e ] kt2= %.3e x= %.3e value= %.3e \tsum=%.3e,%.3e %.3e\n" , imin,imax,pow(par[0],2),par[1], val,sum,total,sum-total);
		imin=imax;
	}
	
	
	Kahn_clear(acc);
	sum=0;
	imin=1.0e-8;
	imax=10*PI/(4*par[0]);
	printf("\n***************************\nexact\n****************************\n");
	for (int i=0;i<10;++i){
		imax+=scale;
		val=dclenshaw<testfunc,double*>(cc, func2, par,imin,imax,1.0e-10,1.0e-10);
		sum+=val;
		acc+=val;
		total=Kahn_total(acc);
		printf("[%.3e %.3e ] kt2= %.3e x= %.3e value= %.3e \tsum=%.3e,%.3e %.3e\n" , imin,imax,pow(par[0],2),par[1], val,sum,total,sum-total);
		imin=imax;
	}
	/*
	Kahn_clear(acc);
	sum=0;
	imin=1.0e-8;
	imax=PI/(4*par[0]);
	xg.init(1.0e-8,1,0.9,pow(imin,-2)*2,1.0,0.2 );
	xg.set_x(par[1]);
	printf("\n***************************\nchebyshev\n****************************\n");
	for (int i=0;i<10;++i){
		imax+=scale;
		val=dclenshaw<testfunc_cheb,double*>(cc, func1, par,imin,imax,1.0e-10,1.0e-10);
		sum+=val;
		acc+=val;
		total=Kahn_total(acc);
		printf("[%.3e %.3e ] kt2= %.3e x= %.3e value= %.3e \tsum=%.3e,%.3e %.3e\n" , imin,imax,pow(par[0],2),par[1], val,sum,total,sum-total);
		imin=imax;
	}
	*/
		
	return 0;
}
	
