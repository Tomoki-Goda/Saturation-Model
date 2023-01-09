#include<iostream>
#include<cuba.h>
#include<cmath>
#include"./control.h"
#include"./control-default.h"
#include"./constants.h"

#include"./clenshaw-curtis.hh"
extern PREC INT_PREC;


#include"./aF.hh"
#include"./sigma.hh"

#ifndef ALPHA_RUN
	#define ALPHA_RUN 0 
#endif
#ifndef MODX
	#define MODX 0
#endif
#ifndef REUSE
	#define REUSE 0 
#endif
#define TEST 0

inline PREC modx(const PREC x, const PREC Q2, const  PREC mf2){
#if MODX==1
	return( (x*(1+4*mf2/Q2)));
#else 
	return( x);
#endif
}


int F2_integrand_T(const int * ndim, const PREC intv[],const int * ncomp,PREC f[], void* p){
	*f=intv[0]*intv[1]*intv[1];
	return 0;
}

void llTest(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const long long int nvec,
  const cubareal epsrel, const cubareal epsabs,
  const int flags,
  const long long int mineval, const long long int maxeval,
  const int key,
  const char *statefile, void *spin,
  int *nregions, long long int *neval, int *fail,
  cubareal integral[], cubareal error[], cubareal prob[]){
  printf("calling\n");
  
  }

double F2_kt(const PREC x,const  PREC Q2,const  PREC mf2,const PREC (& par)[]){
	static int count;
	std::string type="gbw";
#if R_FORMULA==1
	const int ndim=2;
	Sigma sigma[]={Sigma(type,par) ,Sigma(type,par) ,Sigma(type,par) };
	Integrand_r integrands[]={
		Integrand_r(type,modx(x,Q2,MASS_L2),Q2,MASS_L2,sigma[0]) ,
		Integrand_r(type,modx(x,Q2,MASS_C2),Q2,MASS_C2,sigma[1]) ,
		Integrand_r(type,modx(x,Q2,MASS_B2),Q2,MASS_B2,sigma[2])
       	};
	const int key =13;
#else
	const int ndim=3;
	//create once
	//static F2_integrand F2_integrand_A;
	static Gluon gluon[]={Gluon(type,par) ,Gluon(type,par) ,Gluon(type,par) };
	static Integrand_kt integrands[]={
		Integrand_kt(type, x, Q2, MASS_L2, gluon[0]),
	       	Integrand_kt(type, x, Q2, MASS_C2, gluon[1]),
	       	Integrand_kt(type, x, Q2, MASS_B2, gluon[2])
       	};

	gluon[0].set_par(par);
	gluon[1].set_par(par);
	gluon[2].set_par(par);
	integrands[0].set_kinem(x,Q2,MASS_L2);
	integrands[1].set_kinem(x,Q2,MASS_C2);
	integrands[2].set_kinem(x,Q2,MASS_B2);

	const int key =11;
#endif
	
	const long long int mineval=pow(15,ndim), maxeval=1/pow(INT_PREC/10,2);//use llChure if larger than ~1.0e+9
	const long long int nstart=1.0e+2,nincrease=1.0e+2;
	long long int neval=0;
	//const int mineval=pow(10,ndim), maxeval=1/pow(INT_PREC/10,2);//use llChure if larger than ~1.0e+9
	//const int nstart=1.0e+2,nincrease=1.0e+2;
	//int neval=0;

	const int flag= 0+4*0+8*1+16*0+32*0;
	
	int nregions=0,fail=0;
	//cubareal
	PREC integral[3]={0},error[3]={0},prob[3]={0};
	//int spin=0;
	char statefile[100]="";
 	PREC result=0;
 	//printf("%d: Cuhre x=%.2e Q2=%.2e mf2=%.2e\n",count++, x,Q2,mf2);
	llCuhre(ndim, 1,
	//llTest(ndim, 1,
#if R_FORMULA==1
		&F2_integrand_B,
#else
		&F2_integrand_A,
		//&F2_integrand_A1,
#endif  
		(void*)integrands,
		 1, INT_PREC, INT_PREC/10, flag, mineval, maxeval, key, statefile, NULL, &nregions, &neval,  &fail, integral, error, prob
	);
	//printf("Cuhre-End\n");
	//cubawait(&spin);

/*	llVegas(ndim,  1,
#if R_FORMULA==1
		&F2_integrand_B,
#else
		&F2_integrand_A,
		//&F2_integrand_A1,
#endif  
		(void*)integrands,
		 1, INT_PREC, INT_PREC/10, flag,0,
  		mineval,  maxeval,
 		nstart, nincrease, 1000,
		gridno, NULL, NULL,
 		&neval, &fail,
 		integral, error, prob
	);
  */
  
 	result=Q2/(2*PI) *integral[0];
 
 	if(std::isnan(result)+std::isinf(result)!=0){
 		printf("%.3le encountered \n",(double)result);
 	 	return 0;
 	}
 	if(prob[0]>0.1){
 		printf("Q2=%.3e x=%.3e res=%.3e error=%.3e prob=%.3e neval=%lld nregions=%d\n", Q2, x, integral[0],error[0], prob[0], neval,nregions);
 	}
 	//getchar();
 	return ((double)(result));
}



