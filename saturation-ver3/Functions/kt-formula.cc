#include<iostream>
#include<cmath>
#include"/home/tomoki/Numerics/clenshaw-curtis-gauss-legendre.hh"
//#include"/home/tomoki/Numerics/clenshaw-sector.h"
#include"/home/tomoki/Numerics/clenshaw.h"
#include"/home/tomoki/Numerics/dgauss.h"
//#include"/home/tomoki/Numerics/chebyshev.hh"

#include"./control.h"
#include"./control-default.h"
#include"./constants.h"
extern double INT_PREC;
#include"./Photon.hh"
///////////////////////////////////////////////
// Formulae are found in M. A. Kimber Thesis
///////////////////////////////////////////////
#ifndef ALPHA_RUN
	#define ALPHA_RUN 0 
#endif
#ifndef MODX
	#define MODX 0
#endif
#ifndef REUSE
	#define REUSE 0 
#endif


static double alpha(double mu2){
	if(mu2<2*LQCD2){
		mu2=LQCD2;	
	}
	double b0=(11.0*3.0-3.0*2.0)/12.0;
	double val=1.0/(b0*log(mu2/LQCD2));
	//printf("%.5e\n",val);
	//return 0.2;
	return val;
}

static double aF(double x,double k2,double *par){
	if(x>1){
		x=1;
	}
	if(x>1+1.0e-5){
		printf("x/z= %.5e\n",x);
	}
	double lambda=par[1];
	double x0=par[2];

	double Qs2=pow(x0/x,lambda);
	double val=3.0/(4*PI*PI)*k2/Qs2*exp(-k2/Qs2);



	if(std::isnan(val)==1){
		return(0);
	}
	return (par[0]*val) ;
}

static int I_array(double beta, double kappa_t_prime2,double kt2, double Q2, double mf2, double *I){
	double N1=beta*(1-beta)*Q2+mf2;
	double N2=kappa_t_prime2+pow(1-beta,2)*kt2;
	double N3=kappa_t_prime2-pow(1-beta,2)*kt2;
	double N4=kappa_t_prime2+beta*(1-beta)*kt2;
	
	double nsqrt=sqrt(N1*N1+2*N1*N2+N3*N3);
	I[0]=(N1*N2+N3*N3)/pow(nsqrt,3.0);
	I[1]=(N3-(1-2*beta)*N1)/((N1+N4)*nsqrt);	
	I[2]=(N1+N2)/pow(nsqrt,3.0);
	I[3]=(2*(1-beta))/((N1+N4)*nsqrt);	
	return 0;
}


//static 
inline double inv_z(double beta, double kappa_t_prime2,double kt2, double Q2, double mf2){
	return( 1+(kappa_t_prime2+mf2)/(beta*(1-beta)*Q2)+kt2/Q2);
}

static double F2_integrand(double beta, double kappa_t_prime2,double kt2, double Q2, double mf2, double x,double*par ){
	double val=0;
	double I[4];
	I_array(beta, kappa_t_prime2, kt2, Q2, mf2, I);

	val+=(beta*beta+pow(1-beta,2))*(I[0]-I[1]);
	val+=(mf2+4*Q2*beta*beta*pow(1-beta,2) )*(I[2]-I[3]);
	val*=aF(x*inv_z(beta,kappa_t_prime2,kt2,Q2,mf2), kt2 ,par);
#if ALPHA_RUN==1
	val*=alpha(kappa_t_prime2+kt2+mf2+1)/0.2;
	//printf("%.5e\n", val);
#endif
	return(val);	
}
static double F2_integrand_A(double **param){
	double *kinem=param[0];
	double *par=param[1];
	return( F2_integrand(kinem[0],kinem[1],kinem[2],kinem[3],kinem[4],kinem[5],par) );
}


static double kappa_integrand(double *kappa_t_prime2_c,void*param){
	double kappa=*kappa_t_prime2_c;
	/*if(kappa>1){
		printf("we have a problem\n");
	}*/
	//double jac=1;
	double jac=pow( 1-kappa,-2);
	kappa=kappa/(1-kappa);

	double **par=(double**)param;
	par[0][1]=kappa;
	
	return( F2_integrand_A(par) *jac );
}

Clenshaw_Curtis kappa_integrator(16);
//Chebyshev_Gauss kappa_integrator2(32);
static double beta_integrand(double *Beta,void*param){
	double **par=(double**)param;
	par[0][0]=*Beta;
	double beta	=par[0][0];
	double kt2	=par[0][2];
	double Q2	=par[0][3];
	double mf2	=par[0][4];
	double x	=par[0][5];

	double kap_max=beta*(1-beta)*(Q2*((1-x)/x)-kt2)-mf2;
	//double kap_max=(-(mf2*x) + pow(beta,2)*(Q2*(-1 + x) + kt2*x) + beta*(Q2 - kt2*x - Q2*x))/x;
	double kap_min=1.0e-8;
	if(kap_max<kap_min){
		//printf("%.3e %.3e\n",kap_max,kap_min);
		kap_max=kap_min;
		//return 0;
	}
	kap_max=kap_max/(1+kap_max);
	kap_min=kap_min/(1+kap_min);
	double val=0;
	if(kap_min<1.0e-5){
		val+=dgauss(&kappa_integrand,(void*)param,kap_min,1.0e-5,INT_PREC,INT_PREC/10);
	}
	if(kap_max>1.0e-5){
		kappa_integrator.name="kappa";
		kappa_integrator.DIV=4;
		val+=kappa_integrator.integrate(&kappa_integrand,(void*)param,1.0e-5,kap_max,5, INT_PREC,INT_PREC/10);
	}
	//val=dclenshaw(&kappa_integrand,(void*)param,1.0e-10,kap_max, INT_PREC);
	return(val);

}
Clenshaw_Curtis beta_integrator(16);
static double kt_integrand(double* kt2_c,void* param ){
	double **par=(double**)param;
	double kt2=*kt2_c;
	/*if(kt2>1){
		printf("we have a problem\n");
	}*/
	double jac=pow(1-kt2,-2);
	kt2=kt2/(1-kt2);
	//double kt2	=par[0][2];
	double Q2	=par[0][3];
	double mf2	=par[0][4];
	double x	=par[0][5];
	double kapmin=1.0e-8;
	double quadsqrt=1-4*(kapmin+mf2)/(((1-x)/x)*Q2-kt2);
	if(quadsqrt<0){
		return 0;
	}
	double beta_min=0.5-sqrt(quadsqrt)/2;
	double beta_max= 0.5+sqrt(quadsqrt)/2;
	
	
	
	if(beta_min<1.0e-9){
		beta_min=1.0e-9;
	}
	if((1-beta_max)<1.0e-9){
		beta_max=1-1.0e-9;
	}
	//printf("beta %.3e %.3e\n", beta_min,beta_max);
	par[0][2]=kt2;
	//double val=dclenshaw(&beta_integrand,param, 1.0e-10,1-1.0e-10,1.0e-4);
	//double val=beta_integrator.integrate(&beta_integrand,param, 1.0e-10,1-1.0e-10,1.0e-4);
	double val=0;
	if(beta_min<1.0e-5){
		val+=dgauss(&beta_integrand,(void*)param,beta_min,1.0e-5, INT_PREC, INT_PREC/10);
	}
	if(beta_max>1.0e-5){
//		beta_integrator.Reuse=REUSE;
		beta_integrator.name="beta";
		beta_integrator.DIV=4;
		val+=beta_integrator.integrate(&beta_integrand,(void*)param,1.0e-5,beta_max, 5,/*1.0e-4*/INT_PREC, INT_PREC/10);
	}
	return(jac* val/(kt2));
}
Clenshaw_Curtis kt_integrator(16);
double F2_kt(double x,double Q2,double mf2,double *par){
#if MODX==1
	x=x*(1+4*mf2/Q2);
#endif
	double *param[2];	
	param[1]=par;
	double variables[6]={0,0,0,Q2,mf2,x};
	param[0]=variables;
	//double kt_max=8*sqrt(Q2);
	
	double kt_max=-4*(mf2+1.0e-10)+((1-x)/x)*Q2;
	double kt_min=1.0e-2;
	if(kt_max<kt_min){
		return 0;
	}
	kt_max=kt_max/(1+kt_max);
	kt_min=kt_min/(1+kt_min);
	
	//double val=dclenshaw(&kt_integrand,(void*)param,1.0e-10,kt_max, 1.0e-3 );
	double val;
	//kt_integrator.Print=1;
//	kt_integrator.Reuse=REUSE;
	//kt_integrator.Print=1;
	kt_integrator.name="kt";
	val=kt_integrator.integrate(&kt_integrand,(void*)param,kt_min,kt_max, 3, /*1.0e-3*/INT_PREC*10,INT_PREC );
	return(Q2/(2*PI) *  val);
}


////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

double sigma_gbw(double r,double x,double q2, const double * par){
	double sigma_0 =par[0];
	double lambda	=par[1];
	double x_0	=par[2];
	
	if(x_0<0){//to avoid nan since migrad might give negative x0...
		return 0;
	}

	return( sigma_0*(1-exp( - pow(r * Q0, 2) * pow(x_0/x, lambda)/4)) );	
}
double r_integrand(double*R , void *par){
	double **param=(double**)par;
	double r=*R;
	double val=0;
	
	val=psisq_z_int(r,param[0][3],param[0][4] );
	val*=sigma_gbw(r,param[0][5],param[0][3],param[1] );
	val*=2*PI/r;//r^2 comes from photon wave function. just extracted... 2 pi r is angular integration 
	return(val);	
}
Clenshaw_Curtis r_integrator(32);
double F2_r(double x,double Q2,double mf2,double *par){
	//r_integrator.Print=1;
#if MODX==1
	x=x*(1+4*mf2/Q2);
#endif
	double *param[2];	
	param[1]=par;
	double variables[]={0,0,0,Q2,mf2,x};
	param[0]=variables;
	double val;
	r_integrator.name="r";
	val=r_integrator.integrate(&r_integrand,(void*)param,1.0e-6,60,1, INT_PREC*10, INT_PREC*10 );
	//val=dclenshaw(&r_integrand,(void*)param,1.0e-6,60, 1.0e-3 );
	val*=Q2/(4*PI*PI) ;
	//printf("F2= %f\n",val);
	return(val);
}


