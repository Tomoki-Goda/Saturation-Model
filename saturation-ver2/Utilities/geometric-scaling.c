#include<stdio.h>
#include<string.h>
#include<math.h>
#include<stdlib.h>

//#include"control_tmp.h"
#include"control.h"
#include"control-default.h"
#include"constants.h"
 
#include"./plot.c"
//#include"./critical-line.h"
#include"./tmd-gluon.h"

extern double  mod_x(double,double,int);
//extern double generate_points_crit(double x, double **param );

double generate_points_geo(double r, double** par){
	//printf("%.2e %.2e %.2e\n" ,r, *(*(par)), *(*(par)+1) );
	double val = 0;
	double x;
	double qs=par[0][2];
	x=mod_x(**par, *((*par)+1),0);
	val=SIGMA(r/qs/*/0.1973*/,x,*(*(par)+1),*(par+1),*(par+2)) / ( *(*(par+1)) ) ; //0.1973 for GeV<-> fm
	
	return val;
}

int main(int argc, char ** argv){
	int rlen=50;
	double rarr[rlen+1];
	char file_name[500];
	double param[10];
	double sigpar[10];
	double sudpar[10];
	double x, Q2;
	
	double *par[4];

	read_options(argc,argv,param,&x,&Q2, file_name);
	parameter(param,sigpar,sudpar);
	double var[3];
	var[0]=x;
	var[1]=Q2;
	//*par={x,Q2};
	*(par)=var;
	*(par+1)=sigpar;
	*(par+2)=sudpar;	
	float dum;
	char* end;
	
	for(int i=0 ;i<=rlen;i++){
		*(rarr+i)=pow(10,-2+3.5*((double)i)/rlen)/*0.1973*/;
	}
#if (MODEL==1||MODEL==3)	
	approx_xg(sigpar+1);//generate chebyshev coefficients
#endif	
	////////////////////////////////////////////////////////
	double step=(60.0)/(2*n);
	sample_sigma( sample ,  step,  x, Q2, sigpar,  sudpar);	
	var[2]= saturation(step,sudpar,Q2);
	/////////////////////////////////////////////////////
	//var[2]=sqrt(generate_points_crit(x,par ));
	////////////////////////////////////////////////////
	printf("saturation scale=%.3e\n",par[0][2]);
	
	FILE* file=fopen(file_name,"w");
	if(file==NULL){
		printf(" critical line :: file error. %s.\n", file_name);
	}	
	plot(&generate_points_geo,rarr,rlen+1,par,file);
	fclose(file);
	
	return 0;
}	
