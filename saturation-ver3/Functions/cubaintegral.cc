#include<iostream>
#include<cuba.h>
#include<cmath>
#include"./control.h"
#include"./control-default.h"
#include"./constants.h"
extern double INT_PREC;
//#include"./Photon.hh"



#ifndef ALPHA_RUN
	#define ALPHA_RUN 0 
#endif
#ifndef MODX
	#define MODX 0
#endif
#ifndef REUSE
	#define REUSE 0 
#endif

inline double modx(const double x, const double Q2,const  double mf2){
#if MODX==1
	return( (x*(1+4*mf2/Q2)));
#else 
	return( x);
#endif
}

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

static double aF(const double x,const double k2,const double *par){
	double lambda=par[1];
	double x0=par[2];

	double Qs2=pow(x0/x,lambda);
	double val=3.0/(4*PI*PI)*k2/Qs2*exp(-k2/Qs2);



	if(std::isnan(val)==1){
		return(0);
	}
	return (par[0]*val) ;
}

static int I_array(const double beta,const  double kappa_t_prime2,const double kt2, const double Q2, const double mf2, double *I){
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
inline double inv_z(const double beta, const double kappa_t_prime2,const double kt2, const double Q2, const double mf2){
	return( 1+(kappa_t_prime2+mf2)/(beta*(1-beta)*Q2)+kt2/Q2);
}

static double F2_integrand(const double beta, const double kappa_t_prime2,const double kt2, const double Q2,const  double mf2,const  double x, const double*par ){
	if(1-x*inv_z(beta,kappa_t_prime2,kt2,Q2,mf2) <0){
		return 0;		
	}
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

int change_var(double* var,double* jac,const double min, const double max){
	double scale=(max-min);
	double a1=*var;
	double a2=1-a1;
	double den=scale*a2+a1;
	*jac=pow(scale/den,2);
	*var=(min*(scale)*a2+max*a1)/den;
	return 0;
}

double  F2_integrand_A(double kt2, double kappa_t_prime2,double beta, const double x, const double Q2, const double mf2, const double *par){
//#if MODX==1
//	x=x*(1+4*mf2/Q2);
//#endif
	if((1-x)/x *Q2-4*mf2<=0){
		return(0);
	}
	double jac1,jac2,jac3;
	
	const double ktmax=(1-x)/x *Q2-4*mf2;//,kt2=kt201;
	change_var(&kt2, &jac1,0,ktmax);
	
	//double kprimemax=(1-x)/x *Q2/4-mf2;
	const double kprimemax=ktmax/4 ;
	//double kappa_t_prime2=kappa_t_prime201;//just because it is...
	change_var(&kappa_t_prime2, &jac2,0,kprimemax);
	
	double betamin=sqrt(1-4*(mf2/((1-x)/x*Q2) )), betamax ;
	betamax=(1+betamin)/2;
	betamin=(1-betamin)/2;
	change_var(&beta, &jac3,betamin,betamax);
	if(std::isnan(jac1)==1){
		printf("jac1\n");
	}
	if(std::isnan(jac2)==1){
		printf("jac2\n");
	}
	if(std::isnan(jac3)==1){
		printf("Q2= %.3e, x=%.3e, mf2=%.3e sqrt : %.3e jac3\n",Q2,x,mf2,1-4*(mf2/((1-x)/x*Q2) ));
		
	}
	return( jac1*jac2*jac3* F2_integrand(beta,kappa_t_prime2,kt2,Q2,mf2,x,par )/kt2 );
}

int F2_integrand_A(const int *ndim, const double* intv,const int *ncomp,double* f, void* p){
	double **param=(double**)p;
	double beta=intv[0];
	double kappa_t_prime2=intv[1];
	double kt2=intv[2];
		
	const double *kinem=param[0];
	const double *par=param[1];
	
	const double Q2=kinem[0];
	//const double mf2=kinem[1];
	const double x=kinem[2];
	double val=0;

	val+= (2.0/3.0)*F2_integrand_A( kt2, kappa_t_prime2, beta, modx(x,Q2,MASS_L2), Q2,MASS_L2,par);
	val+= (4.0/9.0)*F2_integrand_A( kt2, kappa_t_prime2, beta, modx(x,Q2,MASS_C2), Q2,MASS_C2,par);
	val+= (1.0/9.0)*F2_integrand_A( kt2, kappa_t_prime2, beta, modx(x,Q2,MASS_B2), Q2,MASS_B2,par);
	//std::cout<<*f<<std::endl;
	if(std::isnan(val)+std::isinf(val)!=0){
		printf("evaluation failure %.5e\n", val);
	}
	*f=Q2/(2*PI) *val;
	return 0;
}
/////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////
double psisq_f (const double z,const double r,const double Q2,const double mass2) {
	double	value;
	double	z_bar =  z*z+(1-z)*(1-z);
	double	Qsq_bar =  z*(1-z)*Q2+mass2;
	double	Qsq2 =  sqrt(Qsq_bar)*r;
	//pow(r,2) is to suppress singularity at r=0, it is compensated by the sigma
	if(Qsq2<1.0e-5){//small er approximation
		value =   (z_bar + ( mass2+ pow(2*z*(1-z),2)* Q2 )*pow(r* log(Qsq2),2) );
	}else{
		double	bessel_k0_2 = pow(std::cyl_bessel_k(0,Qsq2),2);
		double	bessel_k1_2 = pow(std::cyl_bessel_k(1,Qsq2),2);
		value = pow(r,2) * (z_bar * Qsq_bar * bessel_k1_2 + ( mass2 + pow(2*z*(1-z),2)* Q2 ) * bessel_k0_2);
	}
	double result=(3*value)/(2*PI*PI);
	if(std::isnan(result)+std::isinf(result)!=0){
 	 	printf("psiisq_f: %.3e encountered \n",result);
 	 	printf("z %.3e  r %.3e Q2 %.3e mass2 %.3e\n",z,r,Q2,mass2);
 	 	getchar();
 	 	return 0;
 	}
	return(result);	
}

double sigma_gbw(const double r,const double x,const double q2, const double * par){
	double sigma_0 =par[0];
	double lambda	=par[1];
	double x_0	=par[2];
	
	if(x_0<0){//to avoid nan since migrad might give negative x0...
		return 0;
	}
	double result= sigma_0*(1-exp( - pow(r * Q0, 2) * pow(x_0/x, lambda)/4)) ;
	if(std::isnan(result)+std::isinf(result)!=0){
 	 	printf("sigma_gbw: %.3e encountered \n",result);
 	 	printf("sigma_0 %.3e  lambda %.3e x_0 %.3e\n",sigma_0,lambda,x_0);
 	 	getchar();
 	 	return 0;
 	}
	return(result);	
}
double F2_integrand_r(const double r,const double z ,const double x,const double Q2,const double mf2,const double*sigpar){
	if((1-x)/x *Q2-4*mf2<0){
		return(0);
	}
	double val=psisq_f (1.0e-6+z/2,r,Q2,mf2);//symmetric in z<->1-z so half. jacobian and 2 cancel
	val*=sigma_gbw(r,x,Q2,sigpar );
	//val*=2*PI/r  ;//r^2 comes from photon wave function. just extracted... 2 pi r is angular integration 
	return(val);	
}



int F2_integrand_B(const int *ndim, const double* intv,const int *ncomp,double* f, void* p){
	double **param=(double**)p;
	double z=intv[0];
	double r=intv[1];
	
	double Q2=param[0][0];
	//double mf2=param[0][1];
	double x=param[0][2];
/*
#if MODX==1
	x=x*(1+4*mf2/Q2);
#endif	
	*f=F2_integrand_r(r,z,x,Q2,mf2,param[1]);
*/
	double jac=0;//=pow((R_MAX-R_MIN)/((R_MAX-R_MIN)*(1-r01)+r01),2);
	change_var(&r,&jac,R_MIN,R_MAX);
	double res=0;
	res+=(2.0/3.0)*F2_integrand_r(r,z,modx(x,Q2,MASS_L2),Q2,MASS_L2,param[1]);
	res+=(4.0/9.0)*F2_integrand_r(r,z,modx(x,Q2,MASS_C2),Q2,MASS_C2,param[1]);
	res+=(1.0/9.0)*F2_integrand_r(r,z,modx(x,Q2,MASS_B2),Q2,MASS_B2,param[1]);

	*f=jac*(Q2)/(2*PI*r) *res;//r^2 comes from photon wave function. just extracted... 2 pi r is angular integration 

	
	//printf("F2= %f\n",val);
	return(0);
}


double F2_kt(const double x,const  double Q2,const  double mf2,const double* par){
#if R_FORMULA==1
	const int ndim=2;
#else
	const int ndim=3;
#endif
	const long int mineval=pow(25,ndim), maxeval=1.0e+7;//use llChure if larger than ~1.0e+9
	const long int nstart=1.0e+2,nincrease=1.0e+2;

	const int flag= 0+4*0+8*1+16*0+32*0;
	long long int neval;
	int nregions,fail;
	cubareal integral[3],error[3],prob[3];
	
	const double *param[2];
	//double par[]={23,3.0e-4,0.3};
	const double kin[]={Q2,mf2,x};
	param[1]=par;
	param[0]=kin;
	const int key =9;
	
	//char statefile[]="cuba.state";
	char statefile[]="";
 	double result;
 	
	
 	llCuhre(ndim, 1,
#if R_FORMULA==1
		&F2_integrand_B,
#else
		&F2_integrand_A,
#endif 
		(void*)param, 1,INT_PREC,INT_PREC/100, flag, mineval,maxeval, key,statefile,NULL, &nregions, &neval,  &fail, integral, error, prob
	);
	
/*	llVegas(ndim, 1,
#if R_FORMULA==1
			&F2_integrand_B,
#else
			&F2_integrand_A,
#endif 
			(void*)param, 1,INT_PREC,INT_PREC/100,flag, 0,	mineval, maxeval, nstart, nincrease, 1000,50, NULL,NULL, &neval,  &fail, integral+1, error+1, prob+1

	);
*/
	if(fail!=0){
		printf("%d",fail );
	}
	

//	 if(prob[0]>0.1){
//	 	//flag=1;
// 	 	std::cout<<"1. "<<integral[0]<<", "<< error[0]<<",  "<<prob[0]<<std::endl;
// 	 	std::cout<<"1. "<<nregions<<", "<< neval<<",  "<<fail<<std::endl;
//		 Vegas(ndim, 1,
//#if R_FORMULA==1
//			&F2_integrand_B,
//#else
//			&F2_integrand_A,
//#endif 
//			(void*)param, 1,INT_PREC,INT_PREC/100,flag, 0,	mineval, maxeval, nstart, nincrease, mineval/10, mineval/20, NULL,NULL, &neval,  &fail, integral+1, error+1, prob+1
//
//		);
 //		if(prob[1]>0.1){
 //			//std::cout<<"2. "<<integral[1]<<", "<< error[1]<<",  "<<prob[1]<<std::endl;
 //			Suave(ndim,1,
//#if R_FORMULA==1	
//				&F2_integrand_B,
//#else
//				&F2_integrand_A,
//#endif 
//				(void*)param, 1, INT_PREC,INT_PREC/100,flag,0, mineval, maxeval,1000,5,25, NULL, NULL,&nregions, &neval, &fail, integral+2, error+2, prob+2
//			);
//			
// 			if(prob[2]>0.1){
 //				std::cout<<"Cuhre. "<<integral[0]<<", "<< error[0]<<",  "<<prob[0]<<std::endl;
 //				std::cout<<"Vegas. "<<integral[1]<<", "<< error[1]<<",  "<<prob[1]<<std::endl;
 //				std::cout<<"Suave. "<<integral[2]<<", "<< error[2]<<",  "<<prob[2]<<std::endl;
// 			}else{
 //				result=integral[2];
// 			}
// 		}else{
 //			result=integral[1];
// 		}
// 	}else{
 		result=integral[0];
 //	}
 	
 	 //return 0;
 	 if(std::isnan(result)+std::isinf(result)!=0){
 	 	printf("%.3e encountered \n",result);
 	 	return 0;
 	 }
 	 return result;
}



