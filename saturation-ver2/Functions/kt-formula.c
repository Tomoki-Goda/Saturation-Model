//#include<iostream>
//#include<fstream>
//#include<cmath>
//#include<minuit.h>
//#include<getopts.h>
//

#include<stdio.h>
#include<math.h>
#include"clenshaw.h"
#include"./control.h"
#include"./control-default.h"
#include"./constants.h"
///////////////////////////////////////////////
// Formulae are found in M. A. Kimber Thesis
///////////////////////////////////////////////
#define ALPHA_RUN 1 

static double alpha(double mu2){
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



	if(isnan(val)==1){
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


static double inv_z(double beta, double kappa_t_prime2,double kt2, double Q2, double mf2){
	double val=1+(kappa_t_prime2+mf2)/(beta*(1-beta)*Q2)+kt2/Q2;
	return val;
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
	double jac=pow( 1-kappa,-2);
	kappa=kappa/(1-kappa);

	double **par=(double**)param;
	par[0][1]=kappa;
	
	return( F2_integrand_A(par) *jac );
}

static double beta_integrand(double *Beta,void*param){
	double **par=(double**)param;
	par[0][0]=*Beta;
	double beta	=par[0][0];
	double kt2	=par[0][2];
	double Q2	=par[0][3];
	double mf2	=par[0][4];
	double x	=par[0][5];

	//double kap_max=beta*(1-beta)*(Q2*((1-x)/x)-kt2)-mf2;
	double kap_max=(-(mf2*x) + pow(beta,2)*(Q2*(-1 + x) + kt2*x) + beta*(Q2 - kt2*x - Q2*x))/x;
	if(kap_max<1.0e-10){
		return 0;
	}
	kap_max=kap_max/(1+kap_max);
	//double kap_min=1.0e-10;
	//kap_min=kap_min/(1+kap_min);
	double val=dclenshaw(&kappa_integrand,param,1.0e-10,kap_max,1.0e-5 );
	return(val);

}
static double kt_integrand(double* kt2_c,void* param ){
	double **par=(double**)param;
	double kt2=*kt2_c;
	/*if(kt2>1){
		printf("we have a problem\n");
	}*/
	double jac=pow(1-kt2,-2);
	kt2=kt2/(1-kt2);
	

	par[0][2]=kt2;
	double val=dclenshaw(&beta_integrand,param, 1.0e-10,1-1.0e-10,1.0e-4);
	return(jac* val/(kt2));
}

double F2_kt(double x,double Q2,double mf2,double *par){
	double *param[2];
	param[1]=par;
	double variables[6]={0,0,0,Q2,mf2,x};
	param[0]=variables;
	double kt_max=8*sqrt(Q2);
	kt_max=kt_max/(1+kt_max);
	double val=dclenshaw(&kt_integrand,(void*)param,1.0e-10,kt_max, 1.0e-3 );
	return(Q2/(2*PI) *  val);
}


/*
int main(int argc, char** argv){
	
	double par[3]={23.3,3.04e-4,0.288};
	//std::fstream file;
	//file.open("f2test.txt",std::fstream::out );
	FILE* file=fopen("f2test.txt","w");

	double x=1.0e-3, mf2=0.0196;
	double val=0;
	double Q2=0;
	//double kappa;	
	//double *params[2];
	//params[1]=par;
	//double param[6]={0.4,0,100,1000,0.0196,1.0e-4};
	//params[0]=param;
	//double kap_max=param[0]*(1-param[0])*(param[3]*((1-param[5])/param[5])-param[2])-param[4];
	double k2;
	for (int i =0; i<21;i++){
		
		Q2=pow(10, -2+5*((double)i)/20);
		val=(2.0/3.0)*F2(x,Q2,0.0196, par) + (4.0/9.0)*F2(x,Q2,0.196,par) + (1/9.0)*F2(x,Q2,21.16,par);
		fprintf(file,"%.5e\t%.5e\n",Q2,val);
		printf("%.5e\t%.5e\n",Q2,val);
		//file<<Q2<<"\t"<<val<<std::endl;
		//std::cout<<Q2<<"\t"<<val<<std::endl;
		
		// k2=pow(10,-2+5*((double)i)/30);
		// val=aF(x,k2,par);
		// file<<k2<<"\t"<<k2*val<<std::endl;
		//kappa=kap_max*((double)i)/50;
		//kappa=kappa/(1+kappa);
		//val=kappa_integrand(&kappa,(void*)params);
		//std::cout<<kappa<<"\t"<<val<<std::endl;
	}
	//file.close();
	fclose(file);

	return 0;


}*/
