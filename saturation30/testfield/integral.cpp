#include<iostream>
#include<cmath>
#include"./simpson1dv1.hh"
#include<ctime>

extern "C" double simps_(double * , double *,
           double *,double *,double* ,double(*)(double*), double *  ,double *,double* ,double *);

double func(double x){
	double p=pow(2.7,-100*pow((x-2),2));
	//double p=pow(x,4);
	return(p);
}
double func2(double* x){
        double p=pow(2.7,-100*pow((*x-2),2));
        return(p);
}


int main(){
	
	static double res=0;
	static double error=0;
	
	static double argmin=0.0;
	static double argmax=4.0;
	static double eps=0.0;
	
	static double dum1,dum2,dum3;
	clock_t t;

	/*simpson int1(6,3,10);
        clock_t t;
        t=clock();
	for(int i=0; i<10;i++){
		int1.sintegrate(&func, argmax ,argmin , eps , eps ,&res,& error);
	}
        t=clock()-t;

        printf("%e seconds\n",( (double)t)/CLOCKS_PER_SEC);
        printf( "Result G %e +/- %.3e\n", res,error);
	*/
	res=0;	
	simpson int2(2,1,5);
	 
	t=clock();
	for(int i=0; i<1;i++){
		int2.sintegrate(&func, argmax ,argmin , eps , eps ,&res,& error);
	}
	t=clock()-t;
	printf("%e seconds\n",( (double)t)/CLOCKS_PER_SEC);
	printf( "Result G %e +/- %.3e\n", res,error);
/*	
	double result=0.0;
	t=clock();
	for(int i=0; i<3;i++){
		simps_(&argmin,&argmax,&eps,&eps,&eps,&func2,&dum1 ,&result ,&dum2,&dum3);
	}
	t=clock()-t;
	printf("%e seconds \n",( (double)t)/CLOCKS_PER_SEC);

	printf("Result S: %e\n",result);	*/
	return 0;
}
