#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include"control.h"
#include"control-default.h"
#include"constants.h"
#include"./plot.c"
//#include"../gluon-chebyshev.h"
//#include"../dipole-cross-section.h"

//#if (MODEL==3||MODEL==22||MODEL==2)
//#include"../sudakov2.h"
//#endif


//Syntax is ./main -in <input dir > -out <output dir > -Q2 <Q2 > -x <x >


extern void approx_xg(double*);
extern double SIGMA(double ,double,double ,double *,double*);
extern double  mod_x(double,double,int);
extern int parameter(double*,double*,double*);



double generate_points(double r, double** par){
	//printf("%.2e %.2e %.2e\n" ,r, *(*(par)), *(*(par)+1) );
	double val = 0;
	double x;
	for(int i =0; i<( NF-1);i++){
		x=mod_x(**par, *((*par)+1),i);
		val+=SIGMA(r/0.1973,x,*(*(par)+1),*(par+1),*(par+2)) / ( *(*(par+1)) ) ; //0.1973 for GeV<-> fm
	}
	return val/(NF-1);
}


int main(int argc , char ** argv){
	char file_name[500];
	int rlen=100;
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
		*(rarr+i)=pow(10,-2+3*((double)i)/rlen)*0.1973;
	}
#if (MODEL==1||MODEL==3)
	approx_xg(sigpar+1);
#endif

	FILE *file=fopen(file_name,"w");
	if(file==NULL){
		printf("DP:: file error. %s.\n", file_name);
	}	
	plot(&generate_points,rarr,rlen+1,par,file);
	fclose(file);
}

