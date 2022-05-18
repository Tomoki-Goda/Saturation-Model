#include<cmath>
#include<iostream>
#include<string>
#include<vector>
#include<fstream>

#include"../main.h"


////////////////////////////////////////////////////////////
////////////////// common functions ////////////////////////
////////////////////////////////////////////////////////////
double mod_x(double x){
	return(x);
}

double alpha_s(double mu2 ){
	double b0= ((double)(33 -2*NF))/(12*PI);
	return( 1/(b0* log(mu2/LQCD2));//LQCD2 lambda_QCD ^2
	
}

//in simps.f
extern "C" double simps_(double *, double *, double *,double *,double* ,double(*)(double*), double *  ,double *,double* ,double *);

/////////////////////////////////////////////////////////////
//////////////////////////// GBW ////////////////////////////
/////////////////////////////////////////////////////////////

double sigma_gbw(double r, double x,double Q2, double * par){
	double sigma_0	=par[0];
	double lambda	=par[1];
	double x_0	=par[2];
	
	double xm=xod(x);

	return( sigma_0*(1-exp( -r*r * pow(x_0/xm, lambda))) );	
}


/////////////////////////////////////////////////////////////
//////////////////////////// BGK ////////////////////////////
/////////////////////////////////////////////////////////////

double sigma_bgk(double r, double x, double Q2, double * par){
	double sigma_0		=par[0];
	double A_g		=par[1];
	double lambda_g	=par[2];
	double C		=par[3];
	double mu02		=par[4];
	
	double xm=xod(x);	
	double mu2=C/(r*r)+mu02;
	double expo = (pow(r * PI,2) * alpha_s(mu2)* xgpdf(xm,mu2))/ (3* sigma_0);
	
	return( sigma_0*(1-exp(-expo))) ;	
}


/////////////////////////////////////////////////////////////
//////////////////////////// GBS ////////////////////////////
/////////////////////////////////////////////////////////////
class gbs_int{
	private:
		double Q2;
		double x;
		double *par;
		double R;
		
		double x_0=par[2];
		double lambda=par[1];

		double Qs2 =pow(x_0/x, lambda);
		
		double integrand(double *r_ptr){
			double r=*r_ptr;
			double laplacian_sigma=-r *log(R/r)*exp(-Qs2*pow(r,2) /4)*Qs2*(1-(Qs2*pow(r,2))/4);
			return sudakov(r,Q2) *laplacian_sigma;
		}
		
	public:
		gbs_int(double x,double Q2,double R,double *par);
	
		double integrate(){
			double res=0;
			double rmi=1.0e-6;
			double rma=1.0e+3;
			double eps=1.0e-6;
			double dum1,dum2,dum3;

			simps_(&rmi, &R, &eps,&eps,&eps, &integrand ,&dum1,&res,&dum2,&dum3 );
				
			return(res);
		}
}; 

double sigma_gbs(double r, double x, double Q2, double * par){
	gbs_int sigma(x,Q2,r,par);
	return(sigma.integrate());
}

///////////////////////////////////////////////////////////
////////////////// all together //////////////////////////
//////////////////////////////////////////////////////////

double sigma(double  r, double x, double Q2, double * par,unsigned model){
	double (*funcptr)( double , double , double , double *);
	switch (moddel){
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
	return( (*funcptr)(r,x,Q2,*par));
	
}

