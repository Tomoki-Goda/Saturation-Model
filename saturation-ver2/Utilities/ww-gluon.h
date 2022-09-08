#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<cfortran.h>

#include"control.h"
#include"control-default.h"
#include"constants.h"

#include"./Parameters.h"

//#define PHI 1

#include"./plot.c"
//#include"./tmd-gluon.h"

extern double KBN_sum(const double *arr,int len);
//#include"./kahnsum.h"
//#include"./critical-line.h"

//#include"./main.h"
extern double dbesj0_(double*);
extern double dbesj1_(double*);
extern double mod_x(double, double,int);
extern double SIGMA(double , double ,double , const double *, const double*);
extern double BASE_SIGMA(double , double ,double , const double *);

extern double dgauss_(double (*)(const double*), double*,double*,double *  );
extern double dgquad_(double (*)(const double*),  const double*, const double*, const int*  );
extern double dadapt_(double(* )(const double*),double*,double*,int*,double*,double* ,double*,double*);

extern int parameter(const double*,double*,double*);
extern void approx_xg(const double *);

///////////////////////////////////////////////////////////
static double X, K,Q2,SIGPAR[10], SUDPAR[10];   //////////
///////////////////////////////////////////////////////////

double ww_integrand(const double * r){
	double R=*r;
	//double jac=1.0/pow(1-R,2);
	//R=R/(1-R);
	double kr=K*R;
	double bes=dbesj0_(&kr);
	double val=2*PI*bes* SIGMA(R,X,Q2,SIGPAR,SUDPAR)/R;
	return val;

}

double ww_integrand_grad(const double * r){
	double R=*r;
	//double jac=1.0/pow(1-R,2);
	//R=R/(1-R);
	//printf("%f %f\n",R,jac);
	double kr=K*R;
	double bes1=dbesj1_(&kr);
	double bes0=dbesj0_(&kr);
	double val=2*PI*(bes0-kr*bes1)* SIGMA(R,X,Q2,SIGPAR,SUDPAR)/R;
	//printf("%.3e\n",val);
	return val;

}

double ww_integral(){
	double res,err;
	double max=(2)*PI*25/K, min=1.0e-5;
	//int n=96;
	//res=dgquad_(&ww_integrand,&min,&max,&n);
	
	double n=1.0e-5;
	res=dgauss_(&ww_integrand,&min,&max,&n);
	
	return res;
}

double ww_grad(){

	double res,err;
	double max=(2)*PI*25/K, min=1.0e-5;
	//int n=96;
	//res=dgquad_(&ww_integrand,&min,&max,&n);
	
	double n=1.0e-7;
	res=dgauss_(&ww_integrand_grad,&min,&max,&n);
	
	return res;
}
///////////////////////////////////////////////////////////////////////////////////////////
double ww_saturation(){
	double k_step=0.1;
	double k_min=0.1;
	//double k=k_min;
	double prev=1;
	double val;
	double step=(60.0)/(100);
	
	for(int i=0;i<6;i++){
		//printf("%.3e\t%.3e\n",k_min,k_step );
		K=k_min;
		
		for(int j=0;j<100;j++){
			val=ww_grad();
			//printf("%d : %.3e, %.3e, %.3e, %.3e\n",j ,K,k_min, val, prev );
			if(j!=0 ){
				if(prev*val<0){
					k_min=K-k_step;
					//printf("%d : %.3e, %.3e, %.3e, %.3e\n",i ,k,k_min, val, prev );
					break;
				}
			}
			if(fabs(val)>fabs(prev)){
				//printf("error\t%.3e\t%.3e\n",val,prev);
			}
			//printf("%f, %f, %f\n",k, val, prev );
			prev=val;
			K+=k_step;
		}
		if(fabs(prev)<1.0e-5){
			break;
		}
		if(i!=5){
			k_step*=0.05;
		}
		
	}
	//printf("\nk=%.3e, val= %.3e, prev= %.3e\n\n",K, val, prev );
	
	K=K - k_step*fabs(val/(prev-val) );//weighted mid point 
	return(K);
}
