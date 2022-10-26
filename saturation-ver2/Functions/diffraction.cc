#include<cmath>
#include<iostream>
#include<fstream>

#include"constants.h"
#include"control.h"
#include"control-default"


//double      cyl_bessel_k( double v, double x );
////////////////////////////////////////////////////////////////////////////
//Formulae can be found in J. R. Forshaw, et al. Nucl. Phys. A675, (2000) //
////////////////////////////////////////////////////////////////////////////
#define PI 3.1415

extern "C" double SIGMA(double ,double, double,double*,double*);
extern "C" int parameter(double*, double*,double*);
extern "C" int approx_xg(double*);
extern "C" double dgquad_(double (*func)(double*),double*, double *, int*);

struct impact{
	int index;
	double z,ep,k,beta,xp,Q2, mf2, kt2;
	double *sigpar,*sudpar;
	//double function(double,double,double,double*);
	//double function(double,double,double,double*,double*);
	void set_extern(double a2, double a3,double a4,double a5){
		//index=a1;
		beta=a2;
		xp=a3;
		Q2=a4;
		mf2=a5;
		}
	
	void set_z(double a1){
		z=a1;
		ep=sqrt(z*(1-z)*Q2+mf2);
		k=sqrt(z*(1-z)*Q2*(1-beta)/beta -mf2);
	}
};
/////////////////////////////////////
//////  stores global variables /////
static struct impact param;
////////////////////////////////////

double phi_integrand(double *R){
	double r=*R;
	double k=param.k, ep=param.ep, xp=param.xp, Q2=param.Q2, index=param.index;

	double val=r*std::cyl_bessel_j(index,std::sqrt(k)*r)*std::cyl_bessel_k(index, std::sqrt(ep)*r);//	*(param.func)(r,x,Q2,sigpar);
	double x=xp;//Q2*(1/xp-1);
	val*=SIGMA(r,x,Q2,param.sigpar,param.sudpar);
	return(val);
}

double phi(int index,double z){
	param.index=index;
	param.set_z(z);

	int N=96;
	double val,min=1.0e-4,max=1.0e+4;
	val=dgquad_(&phi_integrand,&min,&max,&N);
	return(val);
}

double phi2_integrand_u(double *U){
	double u=*U;
	double kt=std::sqrt(param.kt2) , z=param.z, xp=param.xp, Q2=param.Q2;
	int index=2;
	
	double val=u*std::cyl_bessel_k(index,std::sqrt(z/(1-z))*u)*std::cyl_bessel_j(index, u );//	*(param.func)(r,x,Q2,sigpar);
	double x=xp;//Q2*(1/xp-1);
	val*=SIGMA(u/kt,x,Q2,param.sigpar,param.sudpar);
	return(val);
}
double phi2_integrand_kt(double *K){
	double kt2=*K;
	param.kt2=kt2;
	double min=1.0e-4, max=1.0e+4;
	int N=96;
	double phi2=dgquad_(&phi2_integrand_u,&min,&max,&N );
	double z=param.z, Q2=param.Q2;
	double val=phi2*phi2*std::log((1-z)*Q2/(kt2));
	return(val);
}

////////////////////////////////////////////////////////////////////
double FD_L_integrand(double *Z){
	double z=*Z;
	double phi0=phi(0,z);
	double val=phi0*phi0*pow(z * (1-z),3);
	return(val);
}
double FD_T_integrand(double *Z){
	double z=*Z;
	double phi0=phi(0,z);
	double phi1=phi(1,z);
	//!!! NOTE !!!, in the phi function parameters k and ep are modified as they depend on z.
	double val=pow(param.ep,2)*(z*z+(1-z)*(1-z))*phi1*phi1;
	val+=param.mf2*phi0*phi0;
	val*=z*(1-z);
	return(val);
}	
double FD_g_integrand(double *Z){
	double z=*Z;
	//param.index=index;
	param.set_z(z);
	double beta=param.beta;
	int N=96;

	double val,min=1.0e-4,max=(1-z)*param.Q2;
	val=dgquad_(&phi2_integrand_kt ,&min,&max,&N);
	val*=(pow(1-beta/z,2)+pow(beta/z,2) )/pow(1-z,3);
	return(val);
}

////////////////////////////////////////////////////////////////////
double FD_LT(int pol){
	//!!!!!!!!  Run param.set_extern(beta,xp,Q2,mf2) before use !!!!!!!!!!!
	const double hc22=pow(0.3893,2), BD=6.0, alpha_s=0.2;//Unknown hc2 inherited.
	double factor, min,max;
	double (*funcptr)(double*);
	switch(pol){
		case 't':
			factor=2*3*pow(param.Q2,2)/(hc22*128*pow(PI,4)*param.beta*BD);
			funcptr=&FD_T_integrand;
			max=0.5;
			//Half the regon. cf. GBW int over k for Q2(1-beta)/(4*beta)>k2>0.
			// sym. z<->1-z. so factor 2.
			min=(1-std::sqrt(1-4*(param.mf2*param.beta/(param.Q2*(1-param.beta))) ) )/2;
			break;

		case 'l':
			factor=2*3*pow(param.Q2,3)/(hc22*32*pow(PI,4)*param.beta*BD);
			funcptr=&FD_L_integrand;
			max=0.5;
			min=(1-std::sqrt(1-4*(param.mf2*param.beta/(param.Q2*(1-param.beta))) ) )/2;

			break;
		case 'g':
			factor=81*param.beta*alpha_s/(512*pow(PI,5)*hc22*BD);
			funcptr=&FD_g_integrand;
			max=1;
			min=param.beta;
			break;

		default:
			printf("error:: choose polarizaion\n");
		}

	int N=96;
	double res=dgquad_(funcptr,&min,&max, &N);
	//res*=2;//see the comm. above.

	res*=factor;
	res*=2.0/3.0;//sum of charges uds;
	return(res);
}

double FD_T(double beta,double xp,double Q2,double mf2){
	param.set_extern(beta,xp,Q2,mf2);
	double res=FD_LT('t');
	return(res);
}
double FD_L(double beta,double xp,double Q2,double mf2){
	param.set_extern(beta,xp,Q2,mf2);
	double res=FD_LT('l');
	return(res);
}
double FD_g(double beta,double xp,double Q2,double mf2){
	param.set_extern(beta,xp,Q2,mf2);
	double res=FD_LT('g');
	return(res);
}



int main(){
	double beta=0.2;
	double Q2=8.5;
	double mf2=MASS_L2;
	double xp;
	double val;
	std::fstream file;
	file.open("testplot.txt",std::fstream::in);
	char type[3]={'t','l','g'};

	for(int i=0;i<50;i++){
		xp=pow(10,-4+3*((double)i)/50);
		param.set_extern(beta,xp,Q2,mf2);
		val=0;
		for(int j=0;j<3;j++){
			val+=FD_LT(type[j]);
		}
		file<<xp<<"\t"<<val<<std::endl;
	}
	file.close();



}
