#include<cmath>
#include<iostream>
#include<fstream>

//double      cyl_bessel_k( double ν, double x );

extern "C" double SIGMA(double ,double, double,double*,double*);
extern "C" double dgquad(double (*func)(double),double*, double *, int*);

struct impact{
	int index;
	double k,beta,xp,Q2;
	double *sigpar,*sudpar;
	//double function(double,double,double,double*);
	//double function(double,double,double,double*,double*);
	void set(int a1,double a2, double a3,double a4,double a5){
		index=a1;
		k=a2;
		beta=a3;
		xp=a4;
		Q2=a5;
	}
};
/////////////////////////////////////
//////  stores global variables /////
static struct impact param;
////////////////////////////////////

double phi_integrand(double *R){
	double r=*R;
	double k=param.k, x=param.x, Q2=param.Q2, deta=param.beta, index=param.index;

	double val=r*std::cyl_bessel_j(k*r)*std::cyl_bessel_k(std::sqrt(beta/(1-beta))*k*r);//	*(param.func)(r,x,Q2,sigpar);
	val*=SIGMA(r,x,Q2,param.sigpar,param.sudpar);
	return(val);
}

double phi(int index,double k, double beta,double xp,double Q2){
	param.set(index,k,beta,xp,Q2);
	int N=96;
	double val,min=1.0e-4,max=1.0e+4;
	val=dgquad(&phi_integrand,&min,&max,&N);
	return(k*k*val);
}

double FD_L_integrand(double *z){
	double 
}


