#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<cfortran.h>

#include"control.h"
#include"control-default.h"
#include"constants.h"

#include"./Parameters.h"


//#define PHI 0

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
	
	int gridsize[3]={100,100,100};
#if (SUDAKOV==0)
	gridsize[2]=1;
#endif
	for(int i=0;i<gridsize[0]; i++){
		x=pow(10,-7+6*((double)i)/gridsize[0]);
		//x=1.0e-4;
		for(int j=0;j<gridsize[1];j++){
			k=pow(10,-4+8*((double)j)/gridsize[1]);
			for(int l=0;l<gridsize[2];l++){
				Q2=pow(10, -2+5*((double)l)/gridsize[2] );
				//Q2=5;
				//sample_sigma( sample ,  step,  x, Q2, sigpar,  sudpar);
				val=af(x,k,Q2,sigpar,sudpar);
				//printf("%.5e\n",val);
				//val=fill_arr(k, step,sudpar,Q2);
				val*=3.0/(4*PI);
				val/=(2*PI*2*PI);
#if (SUDAKOV==0)	
				fprintf(file,"%.10e\t%.10e\t%.10e\n",log(x),2*log(k), val);
				//fprintf(file,"%.10e\t%.10e\n",k*k, val);
#else
				fprintf(file,"%.10e\t%.10e\t%.10e\t%.10e\n",log(x),2*log(k),log(Q2), val);
#endif
			}
		}
	
	}
	fclose(file);
	
	return 0;
	
}



