#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<cfortran.h>

#include"control.h"
#include"control-default.h"
#include"constants.h"

#include"./Parameters.h"


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
	double step=(60.0)/(2*n);
	
	read_options(argc,argv,param,&x,&Q2, file_name);
	parameter(param,sigpar,sudpar);
	printf("%.3e %.3e\n",x, Q2);
	
#if (MODEL==1||MODEL==3)	
	approx_xg(sigpar+1);//generate chebyshev coefficients
#endif

	
	file=fopen(file_name,"w");

	if(file==NULL){
		printf("tmd-gluon:: file can't be opened. %s\n",file_name);
		return 1;
	}
	for (int i=0; i<=20; i++){
		x= pow(10,-6+((double)4*i)/20);
		sample_sigma( sample ,  step,  x, Q2, sigpar,  sudpar);
		
		val= saturation(step,sudpar,Q2);
		//val*=k*k;
		//printf("%.5e\t%.5e\t%.5e\n",x, val, grad_k(val,step));
		
		fprintf(file,"%.5e\t%.5e\n",x, val*val);
	}
	fclose(file);
	
	return 0;
}

