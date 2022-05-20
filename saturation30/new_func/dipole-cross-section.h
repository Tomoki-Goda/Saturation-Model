/*#include<cmath>
#include<iostream>
#include<string>
#include<vector>
#include<fstream>
*/
#include<stdio.h>
#include<math.h>
//#include"../main.h"
#include"./constants.h"
#include"./intfordp.h"


#define BGK 0 
#define UNIFY 0
////////////////////////////////////////////////////////////
////////////////// common functions ////////////////////////
////////////////////////////////////////////////////////////
double mod_x(double x){
	return(x);
}

double alpha_s(double mu2 ){
	double b0= ((double)(33 -2*NF))/(12*PI);
	return( 1/(b0* log(mu2/LQCD2)));//LQCD2 lambda_QCD ^2
	
}

//in simps.f
//extern "C" double simps_(double *, double *, double *,double *,double* ,double(*)(double*), double *  ,double *,double* ,double *);

/////////////////////////////////////////////////////////////
//////////////////////////// GBW ////////////////////////////
/////////////////////////////////////////////////////////////

double sigma_gbw(double r, double x,double Q2, double * par){
	double sigma_0	=par[0];
	double lambda	=par[1];
	double x_0	=par[2];
	
	double xm=mod_x(x);

	return( sigma_0*(1-exp( -r*r * pow(x_0/xm, lambda))) );	
}


/////////////////////////////////////////////////////////////
//////////////////////////// BGK ////////////////////////////
/////////////////////////////////////////////////////////////
#if BGK
double sigma_bgk(double r, double x, double Q2, double * par){
	double sigma_0		=par[0];
	double A_g		=par[1];
	double lambda_g		=par[2];
	double C		=par[3];
	double mu02		=par[4];
	
	double xm=mod_x(x);	
	double mu2=C/(r*r)+mu02;
	double expo = (pow(r * PI,2) * alpha_s(mu2)* xgpdf(xm,mu2))/ (3* sigma_0);
	
	return( sigma_0*(1-exp(-expo))) ;	
}
//#else
//double sigma_bgk(double r, double x, double Q2, double * par);
#endif

/////////////////////////////////////////////////////////////
//////////////////////////// GBS ////////////////////////////
/////////////////////////////////////////////////////////////
double sudakov(double r, double mu2,double* par) {
	//pertubative sudakov
	double C=*(par);
	double mu02=*(par+1);

	double mub2=C/(pow(r,2)) + mu02;
	if (mu2 < LQCD2 || mub2 < LQCD2|| mu2 < mub2) {
		return(0.0);
	}
	
	double b0 = (11*CA-2*NF)/12;
	double val = CA/(2*b0*PI)*(log(mu2/LQCD2)*log(log(mu2/LQCD2)/log(mub2/LQCD2))-log(mu2/mub2));
    return val;
}


double integrand_gbs(double r, double **par ){
	//for the integral integrand has to be in the form 
	//func(double , double**) where integration is over the first argument. second are constants to be passed.
	//hence the following constant  parameters.
	double R=*(*par);
	double x=*(*par+1);
	double xm=mod_x(x);
	double Q2=*(*par+2);
	//double sigma_0=*(*(par+1) );
	double lambda=*(*(par+1)+1);
	double x_0   =*(*(par+1) +2);
	double *sudpar=(*(par+1)+3);//whatever parameter sudakov takes...

	double Qs2 =pow(x_0/xm, lambda);
	double laplacian_sigma=-r *log(R/r)*exp(-Qs2*pow(r,2) /4)*Qs2*(1-(Qs2*pow(r,2))/4);
            
	return( sudakov(r,Q2,sudpar) *laplacian_sigma);
}


double sigma_gbs(double r, double x, double Q2, double * par){
	double param[3]={r,x,Q2};
	double *param_ptr[2]={param, par};
	double result=0;
	simpson1d(&integrand_gbs, param_ptr,0,r,&result);
	return(result);
}


/*
class gbs_int{

//create an object gbs_int. 
// then obj.sigma_gbs(r,x,Q2,par) return the value.
	private:
		double Q2;
		double x;
		double *par;
		double R;
		
		double x_0=par[2];
		double lambda=par[1];
		double *sudpar=(par+3);//for sudakov 4th and 5th

		double Qs2 =pow(x_0/x, lambda);
		
		double integrand(double *r_ptr){
			double r=*r_ptr;
			double laplacian_sigma=-r *log(R/r)*exp(-Qs2*pow(r,2) /4)*Qs2*(1-(Qs2*pow(r,2))/4);
			return( sudakov(r,Q2,sudpar) *laplacian_sigma);
		}
		
	public:
		gbs_int(double X,double q2,double r, double *param){
				set_variables( X, q2,r,  *param)
		}
		void set_variables(double X,double q2,double r, double *param){
			//set global
			x=X;
			Q2=q2;
			R=r;
			par=param;			
		}
	
		double integrate(){
			double res=0;
			double rmi=1.0e-6;
			double rma=1.0e+3;
			double eps=1.0e-6;
			double dum1,dum2,dum3;

			simps_(&rmi, &R, &eps,&eps,&eps, &integrand ,&dum1,&res,&dum2,&dum3 );
				
			return(res);
		}
		double sigma_gbs(double r, double x, double Q2, double * par,){
			sigma_gbs_i.set_variables(x,Q2,r,par);
			return(sigma_gbs_i.integrate());
		}

}; 
*/

///////////////////////////////////////////////////////////
////////////////// all together //////////////////////////
//////////////////////////////////////////////////////////
#if UNIFY
double sigma(double  r, double x, double Q2, double * par,unsigned model){
	double (*funcptr)( double , double , double , double *);
	switch (model){
		case 0:
			funcptr=&sigma_gbw;
			break;
		case 1:
			funcptr=&sigma_bgk;
			break;
		case 2:
			funcptr=&sigma_gbs;
			break;
		}
	return( (*funcptr)(r,x,Q2,par));
	
}
#endif

