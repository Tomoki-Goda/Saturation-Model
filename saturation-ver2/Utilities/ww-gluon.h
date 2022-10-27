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
static struct ww_parameters{
	double X, K,Q2, *SIGPAR, *SUDPAR;
} ww_par;
///////////////////////////////////////////////////////////
///////////////////in sudakov.h//////////////////
extern double sudakov_p(double, double ,const double*);
extern double sudakov_np(double, double, double ,const double*);
extern int compute_mu2(double,const double*,double*,int);
/////////////////////////////////////////////////
//
double ww_integrand(const double * r){
	double R=*r;
	//double jac=1.0/pow(1-R,2);
	//R=R/(1-R);
	double kr= ww_par.K * R;
	double bes=dbesj0_(&kr);
	double q2=ww_par.Q2;

	double sigma=BASE_SIGMA(R, ww_par.X, ww_par.Q2, ww_par.SIGPAR)/(*(ww_par.SIGPAR));
	double val=2*PI*bes*(2*sigma-sigma*sigma)*(*(ww_par.SIGPAR))/R;

	double sud=1;
#if SUDAKOV>=1
	double  mu2;
	compute_mu2(R,ww_par.SUDPAR,&mu2,1);
	sud*=exp(-sudakov_p(mu2,q2,ww_par.SUDPAR));
#if SUDAKOV==2
	sud*=exp(-sudakov_np(R,mu2,q2,ww_par.SUDPAR));
#endif
#endif
	val*=sud;
	return val;

}

double ww_integrand_grad(const double * r){
	double R=*r;
	printf("DISCONTINUED RESULT IS WRONG!!!!\n");
	//double jac=1.0/pow(1-R,2);
	//R=R/(1-R);
	//printf("%f %f\n",R,jac);
	double kr=ww_par.K*R;
	double bes1=dbesj1_(&kr);
	double bes0=dbesj0_(&kr);
	double val=2*PI*(bes0-kr*bes1)* SIGMA(R,ww_par.X,ww_par.Q2,ww_par.SIGPAR,ww_par.SUDPAR)/R;
	//printf("%.3e\n",val);
	return val;

}

double ww_integral(){
	double res,err;
	//double max=(2)*PI*150/ww_par.K, min=1.0e-5;
	double max=100, min=1.0e-5;
	//int n=96;
	//res=dgquad_(&ww_integrand,&min,&max,&n);
	
	double n=1.0e-8;
	res=dgauss_(&ww_integrand,&min,&max,&n);
	printf("%.5e\n",res);	
	return res;
}

double ww_grad(){

	double res,err;
	double max=(2)*PI*25/ww_par.K, min=1.0e-5;
	//int n=96;
	//res=dgquad_(&ww_integrand,&min,&max,&n);
	
	double n=1.0e-7;
	res=dgauss_(&ww_integrand_grad,&min,&max,&n);
	
	return res;
}
///////////////////////////////////////////////////////////////////////////////////////////
/*double ww_saturation(){
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
*/
