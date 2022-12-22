#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include"control.h"
#include"control-default.h"
#include"constants.h"
//#include"../gluon-chebyshev.h"
//#include"../dipole-cross-section.h"
//#include"../photon-wave-function-2.h"
//#include"../simpson-integral.h"
#include"./plot.c"


extern void approx_xg(const double *);
extern int parameter(const double*,double *, double*);
//extern void simpson1dA(double(*)(double, double**),double**,double,double,int,double*,double*); 
extern double SIGMA(double , double ,double ,double *,double*);
extern double psisq_z_int(double, double ,int);
extern double mod_x(double,double, int);
#include"./f2.h"


double f2_slope(double Q2, double**pars){
	double  err;
	*(*(pars)+1)=Q2;
	//printf("x== %.2e\t Q2== %.2e\n", pars[0][0],pars[0][1]);
	//getchar();
	//double xdiff=0.01;
	double xdiff=(**pars)/10; //this is just a choice. choosing the value of 10x, note its \Delta log10(x) so very small.
	double logx=log10(**pars);
	double f2[2];
	//simpson1dA(&f2_integrand,pars,1.0e-5,0.97,500,f2,&err);//f2 is not over r bur r=R/(1-R)
	f2[0]=dclenshaw(&f2_integrand,(void*)pars,1.0e-5,0.97,DGAUSS_PREC);
	(*(*pars)) = pow(10.0, logx+xdiff);
	//simpson1dA(&f2_integrand,pars,1.0e-5,0.97,500,f2+1,&err);
	f2[1]=dclenshaw(&f2_integrand,(void*)pars,1.0e-5,0.97,DGAUSS_PREC);
	//double slope = (log(f2[0]/ f2[1]) )/(log(10)*xdiff);
	double slope = (log10(f2[0]/ f2[1]) )/(xdiff);
	printf("%f\t%f\t%f\n",slope,*(f2),*(f2+1));
	(*(*pars)) = pow(10.0, logx);
	return slope;
}



int main (int argc, char** argv){
	char file_name[500];
	int Q2len=30;
	double Q2arr[Q2len+1];
	double param[10];
	double x=0;
	double Q2=0;
	
	read_options(argc,argv,param,&x,&Q2,file_name );
	//printf("x= %.2e\tQ2= %.2e\n",x,Q2);
	
	for(int i =0;i<=Q2len;i++){
		*(Q2arr+i)=pow(10.0, -2+5*((double)(i))/Q2len);
	}

	double *par[4];
	double var[2];
	var[0]=x;
	var[1]=0;
	
	//printf("x= %.2e\tQ2= %.2e\n",var[0], var[1]);
	//getchar();
	par[0]=var;
	double sigpar[10];
	double sudpar[10];
	parameter(param,sigpar,sudpar);
	par[1]=sigpar;
	par[2]=sudpar;
	double flavour=0;
	par[3]=&flavour;

#if (MODEL==1||MODEL==3)	
	approx_xg(sigpar+1);//generate chebyshev coefficients
#endif

	FILE *file=fopen(file_name,"w");
	if(file==NULL){
		printf("F2-slope:: file error. %s.\n", file_name);
	}	
	plot(&f2_slope,Q2arr,Q2len+1,  par,  file);
	fclose(file);

}



