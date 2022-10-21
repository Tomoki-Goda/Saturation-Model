
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#include"control.h"
#include"control-default.h"
#include"constants.h"

#include"./plot.c"

#include<string.h>



int N_SIMPS=N_SIMPS_R;
//int N_CHEB=N_CHEB_R;

//static double PSI[5][2*N_SIMPS_R+1][MAXN];//pre-evaluated sets of psi
//static double SAMPLES[5][2*N_SIMPS_R+1][MAXN];//sampled points of integrand 
#define MAXN 600
//#include"./kahnsum.h"

extern double KBN_sum(const double *arr,int len);
extern double kahn_sum(const double *arr,int len);

extern double SIGMA_PREC;

extern double SIGMA(double , double ,double , const double *, const double*);
extern double psisq_z_int(double, double ,int);
extern double mod_x(double,double, int);

extern void approx_xg(const double *);

extern int parameter(const double*,double*,double*);


double integrand_sample(double r, double x, double Q2,double* sigpar , double * sudpar){
	double val=0;
	//printf("%f\n",val);
	double xm;
	for(int fl=0;fl<(NF-1);fl++){
		xm=mod_x(x, Q2,fl );
		val+=psisq_z_int(r, Q2, fl) *SIGMA(r,xm,Q2, sigpar,sudpar) ;
	}
	
	return val/(r);//again r everywhere!! read comment in "read-and-fit.h"
}
double generate_points(double r,double **par){
	double x=**par;
	double Q2=*(*par+1);
	double *sigpar=*(par+1);
	double *sudpar=*(par+2);
	double val;
	val=integrand_sample( r,  x, Q2, sigpar, sudpar);
	//printf("%.5e",val);
	return val;
}


int main(int argc , char ** argv){
	char file_name[500];
	int rlen=50;
	double rarr[rlen+1];
	double x, Q2;
	double param[10];
	double sigpar[10];
	double sudpar[10];	
	
	read_options(argc,argv,param,&x,&Q2, file_name);
	parameter(param,sigpar,sudpar);
	double *par[3];
	double var[2];
	*var=x;
	*(var+1)=Q2;
	*(par)=var;
	*(par+1)=sigpar;
	*(par+2)=sudpar;

	for(int i=0 ;i<=rlen;i++){
		//*(rarr+i)=pow(10,-2+3.5*((double)i)/rlen);//*0.1973;
		*(rarr+i)=pow(10,-5+6.5*((double)i)/rlen);//*0.1973;
	}
#if (MODEL==1||MODEL==3)
	approx_xg(sigpar+1);
#endif

	FILE *file=fopen(file_name,"w");
	if(file==NULL){
		printf("Integrand:: file error. %s.\n", file_name);
	}	
	plot(&generate_points,rarr,rlen+1,par,file);
	fclose(file);
}



	


























