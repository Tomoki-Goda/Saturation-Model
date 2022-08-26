#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<cfortran.h>

#include"control.h"
#include"control-default.h"
#include"constants.h"

#include"./Parameters.h"


#include"./plot.c"
#include"./tmd-gluon.h"

int main (int argc, char** argv){

	double val;
	FILE *file;
	char file_name[500];
	double k2, x , Q2;
	double param[10];
	double sudpar[10];
	double sigpar[10];
	double step=((double)R_MAX)/(2*n);
	
	
	
	read_options(argc,argv,param,&k2,&Q2, file_name);
	parameter(param,sigpar,sudpar);
	printf("%.3e %.3e\n",k2, Q2);
		
#if (MODEL==1||MODEL==3)	
	approx_xg(sigpar+1);//generate chebyshev coefficients
#endif
	file=fopen(file_name,"w");
	for(int i=0;i<500; i++){
		x=pow(10,-6+4*((double)i)/500);
		sample_sigma( sample ,  step,  x, Q2, sigpar,  sudpar);
		if(file==NULL){
			printf("tmd-gluon:: file can't be opened. %s\n",file_name);
			return 1;
		}
		val=fill_arr(k2, step,sudpar,Q2);
		//val*=k*k;
		val*=3.0/(4*PI);
			
		fprintf(file,"%.5e\t%.5e\n",x, val);
	}
	fclose(file);
	
	return 0;
	
}



