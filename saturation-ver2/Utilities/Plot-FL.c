

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
//#include"./plot.c"
#include"./options.h"

extern int F_L;
extern void approx_xg(const double *);
extern int parameter(const double*,double *, double*);
extern void simpson1dA(double(*)(double, double**),double**,double,double,int,double*,double*); 
extern double SIGMA(double , double ,double ,double *,double*);
extern double psisq_z_int(double, double ,int);
extern double mod_x(double,double, int);
#include"./f2.h"



int main (int argc, char** argv){
	read_options(argc,argv,&OPTIONS);
	double param[10];
	F_L=1;
	if(((OPTIONS.W)!=((int)0))){
		FIX_W=1;//in f2.h
		//extern, in photon wave function
	}
	FILE * input=fopen(OPTIONS.input_file_name,"r");
	read_parameters(input,param);
	fclose(input);
	double sigpar[10];
	double sudpar[10];
	parameter(param,sigpar,sudpar);
//////////////////////////////////
	double *par[4];
	double var[2];
	par[0]=var;
	par[1]=sigpar;
	par[2]=sudpar;
	
	double flavour=0;
	par[3]=&flavour;

	double val;
	FILE* out=fopen(OPTIONS.output_file_name,"w");
#if (MODEL==1||MODEL==3)	
	approx_xg(sigpar+1);//generate chebyshev coefficients
#endif
	if(FIX_W==1){
		int Q2len=20;
		for(int i =0;i<=Q2len;i++){
			var[1]=pow(10.0, 0+3*((double)(i))/Q2len);
			var[0]=var[1]/(var[1]+pow(OPTIONS.W,2));
			printf("W= %.2e\n",var[0]);
			val=f2(var[1],par);
			fprintf(out,"%.5e\t%.5e\n",var[1],val);		
		}
		
	}else{

		FILE* data=fopen(OPTIONS.data_file_name,"r");
		fscanf(data, "%*[^\n]");
		double x,Q2,f_l,dummy,error;

		while(!feof(data)){
			fscanf(data,"%lf %lf %lf %lf %lf %lf %lf\n",var+1,var,&f_l,&dummy,&dummy,&dummy,&error );
			val=f2(var[1],par);
			fprintf(out,"%.5e\t%.5e\n",var[1],val);
		}
		fclose(data);

	}
	fclose(out);
///////////////////////////////////

}



