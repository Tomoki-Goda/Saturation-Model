#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<cfortran.h>

#include"control.h"
#include"control-default.h"
#include"constants.h"

#include"./Parameters.h"

#define PHI 0 
#include"./plot.c"
#include"./tmd-gluon-2.h"

int main (int argc, char** argv){

	double val;
	FILE *file;
	char file_name[500];
	double k, x , Q2;
	double param[10];
	double sudpar[10];
	double sigpar[10];
	double step=((double)R_MAX)/(2*n);
	
	read_options(argc,argv,param,&x,&Q2, file_name);
	parameter(param,sigpar,sudpar);
	printf("%.3e %.3e\n",x, Q2);
	
#if (MODEL==1||MODEL==3)	
	approx_xg(sigpar+1);//generate chebyshev coefficients
#endif

	//sample_sigma( sample ,  step,  x, Q2, sigpar,  sudpar);
	file=fopen(file_name,"w");

	if(file==NULL){
		printf("tmd-gluon:: file can't be opened. %s\n",file_name);
		return 1;
	}
	for (int i=0; i<100; i++){
		k= pow(10,-1+((double)2*i)/100);
/*#if PHI==1
		val=fill_arr_2(k, step);
#else
		val=fill_arr(k, step,sudpar,Q2);
#endif*/
		val=af(x,k,Q2,sigpar,sudpar);
		//val*=k*k;
		val*=3.0/(4*PI);
		val/=(2*PI*2*PI);
		
		printf("%.5e\t%.5e\n",k*k, val);
		fprintf(file,"%.5e\t%.5e\n",k*k, val);
	}
	fclose(file);
	
	return 0;
}



