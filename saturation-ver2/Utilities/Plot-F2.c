

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
extern void simpson1dA(double(*)(double, double**),double**,double,double,int,double*,double*); 
extern double SIGMA(double , double ,double ,double *,double*);
extern double psisq_z_int(double, double ,int);
extern double mod_x(double,double, int);
#include"./f2.h"


int main (int argc, char** argv){
	char file_name[500];
	int Q2len=30;
	double Q2arr[Q2len+1];
	double param[10];
	double x=0;
	double Q2=0;
	
	read_options(argc,argv,param,&x,&Q2,file_name );
	printf("x= %.2e\tQ2= %.2e\n",x,Q2);

	for(int i =0;i<=Q2len;i++){
		*(Q2arr+i)=pow(10.0, -2+5*((double)(i))/Q2len);
	}

	double *par[4];
	double var[2];
	var[0]=x;
	var[1]=0;
	
	//printf("x= %.2e\tQ2= %.2e\n",var[0], var[1]);
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
		printf("F2:: file error. %s.\n", file_name);
	}	
	plot(&f2,Q2arr,Q2len+1,  par,  file);
	fclose(file);

}



