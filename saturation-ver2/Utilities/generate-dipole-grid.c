#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include"control.h"
#include"control-default.h"
#include"constants.h"
//#include"./plot.c"
#include"options.h"
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
	val=SIGMA(r,x,*(*(par)+1),*(par+1),*(par+2));
	return val;
}


int main(int argc , char ** argv){
	double r, x, Q2, val;
	double param[10];
	double sigpar[10];
	double sudpar[10];	
	
	read_options(argc,argv,&OPTIONS);
	FILE* input_file=fopen(OPTIONS.input_file_name,"r");
	read_parameters(input_file,param);
	fclose(input_file);
	parameter(param,sigpar,sudpar);

#if (MODEL==1||MODEL==3)
	approx_xg(sigpar+1);
#endif

	FILE *file=fopen(OPTIONS.output_file_name,"w");
	if(file==NULL){
		printf("DP:: file error. %s.\n", OPTIONS.output_file_name);
	}	
	int gridsize_x=100;
#if SUDAKOV>=1
	int gridsize_Q2=100;
#else 
	int gridsize_Q2=1;
#endif
	int gridsize_r=100;
	for(int i=0;i<gridsize_x;i++){
		x=pow(10,-7 + 7 * ((double)i)/gridsize_x);
		for(int j=0;j<gridsize_Q2;j++){
			Q2=pow(10,-1+3*((double)j)/gridsize_Q2);
			for(int k=0;k<gridsize_r;k++){
				r=pow(10,-6+ 8*((double)k)/gridsize_r);
				val=SIGMA(r,x,Q2,sigpar,sudpar);
				fprintf(file, "%.5e\t%.5e\t%.5e\t%.5e\n",r,x,Q2,val);
			}
		}
	}

	fclose(file);
}

