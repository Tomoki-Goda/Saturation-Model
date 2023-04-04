#define MODEL 0
#include<cmath>
#include "./Functions/complex.hh"
//#include "cfortran.h"
#include"./Functions/control-default.h"
#include"./Functions/constants.h"
//#include "./Functions/clenshaw.hh"
#include<iostream>
#include<pthread.h>

//#include"./Functions/gluons.hh"


//Collinear_Gluon xgpdf;
extern "C" double xgpdf_(const double* x, const double* QQ,const double* A_g, const double* lambda_g );
extern "C" double xgpdf_integrand_(const double* x, const double* QQ) ;
void* compute(void* arg){
	double *args=(double*)arg;
	double j=1;
	
	//printf("%.3e %.3e %.3e %.3e \n",args[0],args[1],args[2],args[3] );
	//getchar();
	//args[4]=xgpdf(args[0],args[1],args[2],args[3] );
	args[4]=xgpdf_(args,args+1,args+2,args+3);
	//for(int i=0;i<30;i++){
		//j=i;
	//	args[4]=xgpdf_integrand_(&j,args) ;
	//	(args[1]>2)?(printf("2\t\t")):(printf("1\t"));
	//	printf("%.3e \n",args[4]);
	//}
	return NULL;	
}

int main(){
	//pthread_t thread[3];
	double args[5*200];
	int i1;
#pragma opm parallel
{
#pragma omp for
	for(int j=0;j<200;j++){
		args[5*j]=0.0001;//*(j+1);
		args[5*j+1]=pow(10,2);
		args[5*j+2]=1;
		args[5*j+3]=0.01;
		
		//if(j>1) {pthread_join(thread[j-1],NULL); }
		compute((void*)(args+5*j) );
		printf(" %.3e \n", args[5*j+4]);
		
		
	}
	printf(" loop end\n");
}
	
/*	for(int j=0;j<2;j++){
		args[5*j]=0.0001;//*(j+1);
		args[5*j+1]=pow(10,j);
		args[5*j+2]=1;
		args[5*j+3]=0.01;
		
		//if(j>1) {pthread_join(thread[j-1],NULL); }
		i1=pthread_create(thread+j,NULL,compute,(void*)(args+5*j) );
		
	}
*/	for(int j=0;j<2;j++){
		//pthread_join(thread[j],NULL);
		printf(" %.3e \n", args[5*j+4]);
	}
	return 0;
}
