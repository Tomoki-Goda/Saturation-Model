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
#include"./plot.c"
extern void simpson1dA(double(*)(double, double**),double**,double,double,int,double*,double*); 
extern double SIGMA(double , double ,double ,double *,double*);
extern double psisq_z_int(double, double ,int);
extern double mod_x(double,double, int);

extern void approx_xg(double *);
extern int parameter(double*,double *, double*);

double f2_integrand(double r, double ** par){
	double xm;
	double x = *(*par);
	double Q2= *(*(par)+1);
	double *sigpar=*(par+1);
	double *sudpar=*(par+2);
	int fl=(int)(**( par+3)+0.1);

	//parameter(*(par+1),sigpar,sudpar);
	//printf("%f\t%f\t%f\n",par[1][0],par[1][1],par[1][2]);

	double value=0.0;

//	for(unsigned fl=0;fl<1/* NF-1*/;fl++){
//	for(unsigned fl=0;fl< 1;fl++){
		//printf("%d",fl);
	xm=mod_x(x,Q2,fl);
	value+=psisq_z_int(r, Q2, fl)* SIGMA(r,xm,Q2, sigpar,sudpar)/r ;
		
//	}
	//printf("x: %.2e\tx_mod %.2e\t Q2: %.2e\t%f\t%f\t%f\tresult :%f \n", x,mod_x(x,Q2,3),Q2,*(sigpar),*(sigpar+1),*(sigpar+2),value);
	return(value);
}
