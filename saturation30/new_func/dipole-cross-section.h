/*#include<cmath>
#include<iostream>
#include<string>
#include<vector>
#include<fstream>
*/
//#include<stdio.h>
//#include<math.h>
//#include"../main.h"
//#include"./constants.h"
//#include"./simpson-integral.h"
//#include"../gluons.h"

#define BGK 0
#define UNIFY 1 

extern  double xgpdf(double, double);

////////////////////////////////////////////////////////////
////////////////// common functions ////////////////////////
////////////////////////////////////////////////////////////

double alpha_s(double mu2 ){
	double b0= ((double)(33 -2*NF))/(12*PI);
	return( 1/(b0* log(mu2/LQCD2)));//LQCD2 lambda_QCD ^2
	
}//in simps.f
//extern "C" double simps_(double *, double *, double *,double *,double* ,double(*)(double*), double *  ,double *,double* ,double *);

/////////////////////////////////////////////////////////////
//////////////////////////// GBW ////////////////////////////
/////////////////////////////////////////////////////////////

double sigma_gbw(double r, double x,double Q2, double * par){
	double sigma_0	=par[0];
	double lambda	=par[1];
	double x_0	=par[2];
	
	//double xm=mod_x(x);

	return( sigma_0*(1-exp( -r*r * pow(x_0/x, lambda))) );	
}


/////////////////////////////////////////////////////////////
//////////////////////////// BGK ////////////////////////////
/////////////////////////////////////////////////////////////
#if BGK

//extern "C" double xgpdf(double, double);

double sigma_bgk(double r, double x, double Q2, double * par){
	double sigma_0		=par[0];
	double A_g		=par[1];
	double lambda_g		=par[2];
	double C		=par[3];
	double mu02		=par[4];
	
	//double xm=mod_x(x);	
	double mu2=C/(r*r)+mu02;
	double expo = (pow(r * PI,2) * alpha_s(mu2)* xgpdf(x,mu2))/ (3* sigma_0);
	
	return( sigma_0*(1-exp(-expo))) ;	
}
#else
double sigma_bgk(double r, double x, double Q2, double * par){
	return(0.0);
};
#endif

/////////////////////////////////////////////////////////////
//////////////////////////// GBS ////////////////////////////
/////////////////////////////////////////////////////////////
double sudakov(double r, double mu2,double* par) {
	//pertubative+non-perturbative sudakov
	double C= (*(par));
	double r_max=( *(par+1));
	double g1=(*(par+2));
	double g2=(*(par+3));

	double mub2=C*( 1/(pow(r,2)) + 1/pow(r_max,2) ) ;
	if (mu2 < LQCD2 || mub2 < LQCD2|| mu2 < mub2) {
		return(0.0);
	}
	
	double b0 = (11*CA-2*NF)/12;
	double val = CA/(2*b0*PI)*(log(mu2/LQCD2)*log(log(mu2/LQCD2)/log(mub2/LQCD2))-log(mu2/mub2));
	val+=(g1/2) * pow(r,2) + (g2/4) * log(1+pow(r/r_max,2) )*log(pow( mu2 /Q0 ,2) );
    return val;
}


double integrand_gbs(double r, double *par[2] ){
	//for the integral integrand has to be in the form 
	//func(double , double**) where integration is over the first argument. second are constants to be passed.
	//hence the following constant  parameters.
	if(fabs(r)<1.0e-15){
		return(0.0);
	}
	double R=( *(*par) ) ;
	double x=( *(*par+1) );
	//double xm=mod_x(x);
	double Q2=( *(*par+2) ) ;
	//double sigma_0=*(*(par+1) );
	double lambda=( *(*(par+1)+1) );
	double x_0   =( *(*(par+1) +2) );
	double *sudpar=( *(par+1)+3 );//whatever parameter sudakov takes...

	double Qs2 =pow(x_0/x, lambda);
	double laplacian_sigma=r *log(R/r)*exp(-Qs2*pow(r,2) /4)*Qs2*(1-(Qs2*pow(r,2))/4);
       double val=exp(-sudakov(r,Q2,sudpar)) *laplacian_sigma;
	//printf("integrand return for r=%f, %f. \n",r,val);       
	return(val); //sudakov(r,Q2,sudpar) *laplacian_sigma);
}


double sigma_gbs(double r, double x, double Q2, double * par){
	double param[3]={r,x,Q2};

	double *param_ptr[2]={param, par};
	double result=0;
	simpson1d(&integrand_gbs, param_ptr,0.0,r,&result);
	//printf("\n\n");
	return(result);
}


///////////////////////////////////////////////////////////
////////////////// all together //////////////////////////
//////////////////////////////////////////////////////////

//#if UNIFY
/*
 * double sigma(double  r, double x, double Q2, double * par,unsigned model,unsigned char flavour){

	double (*funcptr)( double , double , double , double *);
	double x_mod =mod_x(x,Q2,flavour);

//	switch (model){
//		case 0:
//			funcptr=&sigma_gbw;
//			break;
//		case 1:
//			funcptr=&sigma_bgk;
//			break;
//		case 2:
//			funcptr=&sigma_gbs;
//			break;
//		}

	#if MODEL==0
	funcptr=&sigma_gbw;
	#elif MODEL==1
	funcptr=&sigma_bgk;
	#elif MODEL==2
	funcptr=&sigma_gbs;
	#endif

	return( (*funcptr)(r,x_mod,Q2,par));
	
}
//#endif
*/
