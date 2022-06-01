#include<math.h>
#include<stdio.h>
#include"./simpson-integral.h"



double consts[1]={0.5};

double integrand(double r, double ** param){
	double x=(*(*param));
	double y=(*(*param)+1);
	double a=(*(*(param+1)));

	return(pow(x,a )*exp(-y*pow(r,2)) );
}
double integral(double x,double y ,double *par){
	double var[2]={x,y};
	double *param[2]={var,par};

	double res=0;
	simpson1dA(&integrand, param,0,1,&res);
	return(res);
}

double integral2(double x,double y ,double *par){
	double var[2]={x,y};
	double *param[2]={var,par};

	double res=0;
	simpson1d(&integrand, param,0,1,&res);
	return(res);
}

int main(){
	double res;
	res=integral(2.0,3.0,consts);
	printf("%f\n", res);
	 res=integral2(2.0,3.0,consts);
	printf("%f\n", res);

	return(0);

}

