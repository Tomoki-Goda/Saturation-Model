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
	
	
	
	read_options(argc,argv,param,&k,&Q2, file_name);
	parameter(param,sigpar,sudpar);
	printf("%.3e %.3e\n",k, Q2);
		
#if (MODEL==1||MODEL==3)	
	approx_xg(sigpar+1);//generate chebyshev coefficients
#endif
	file=fopen(file_name,"w");
	if(file==NULL){
			printf("tmd-gluon:: file can't be opened. %s\n",file_name);
			return 1;
	}
	
	for(int i=0;i<50; i++){
		x=pow(10,-7+6*((double)i)/50);
		//sample_sigma( sample ,  step,  x, Q2, sigpar,  sudpar);
/*#if PHI==1		
		val=fill_arr_2(k2, step);
#else
		val=fill_arr(k2, step,sudpar,Q2);
#endif*/

		val=af(x,k,Q2,sigpar,sudpar);
		//val*=k*k;
		val*=3.0/(4*PI);
		val/=(2*PI*2*PI);
			
		fprintf(file,"%.5e\t%.5e\n",x, val);
	}
	fclose(file);
	
	return 0;
	
}



