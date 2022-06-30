#include<stdio.h>
#include<string.h>
#include<math.h>
#include<stdlib.h>

//#include"control_tmp.h"
#include"control.h"
#include"control-default.h"
#include"constants.h"
 
#include"./plot.c"
#include"./critical-line.h"




int main(int argc, char ** argv){
	int xlen=100;
	double xarr[xlen+1];
	char file_name[500];
	double param[10];
	double sigpar[10];
	double sudpar[10];
	double x, Q2;
	
	double *par[3];

	read_options(argc,argv,param,&x,&Q2, file_name);
	parameter(param,sigpar,sudpar);
	double var[2];
	var[1]=Q2;
	//*par={x,Q2};
	*(par)=var;
	*(par+1)=sigpar;
	*(par+2)=sudpar;	
	float dum;
	char* end;
	
#if (MODEL==1||MODEL==3)	
	approx_xg(sigpar+1);//generate chebyshev coefficients
#endif
	//printf("gluon ready\n");
	double point=0;
	for (int i=0;i<=xlen;i++){
		*(xarr+i)=pow(10.0,-6+4*(((double)i)/xlen) );
		
	}
	FILE* file=fopen(file_name,"w");
	if(file==NULL){
		printf(" critical line :: file error. %s.\n", file_name);
	}	
	plot(&generate_points,xarr,xlen+1,par,file);
	fclose(file);
	
	return 0;
}



