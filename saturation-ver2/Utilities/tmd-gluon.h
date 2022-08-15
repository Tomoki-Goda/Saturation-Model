#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<cfortran.h>

#include"control.h"
#include"control-default.h"
#include"constants.h"

#include"./Parameters.h"

#include"./plot.c"
extern double KBN_sum(const double *arr,int len);
//#include"./kahnsum.h"
//#include"./critical-line.h"

//#include"./main.h"
extern double dbesj0_(double*);
extern double dbesj1_(double*);
extern double mod_x(double, double,int);
extern double SIGMA(double , double ,double , const double *, const double*);

extern int parameter(const double*,double*,double*);
extern void approx_xg(const double *);

static const int n=1000;
static double sample[2*1000+1 +2]={0};
static double sample_lap[2*1000+1]={0};
static double summand[2*1000+1]={0};

double sample_sigma(double * sample , double step, double x,double Q2,const double * sigpar, const double * sudpar){ 
	// one only needs to compute a set of sigma once.
	double r=0;
	double xm=0;
	double val=0;
	for(int j=0;j<(2*n+1)+2;j++){
		r=R_MIN+j*step;
		
		val=0;
		//for(int i=0; i<(NF-1); i++){
		for(int i=0; i<1; i++){
			xm=mod_x(x,Q2,i);
			val+=SIGMA(r,xm,Q2,sigpar,sudpar);
			//printf("val=%.3e %.3e %.3e %.3e\n",val,r ,xm,Q2);
		}
		//printf("%.3e\n",val);
		sample[j]=val;
	}
}

double simps_sum(double * sample, int len ,  double step){
	double summand[len];
	double term, val;
	for(int j=0;j<(2*n+1);j++){	
			term=sample[j];
				
			if((j==0)||(j==2*n)){
			
			} else if( (j/2)*2==j ){
				term*=2;
			}
			else{
				term*=4;	
			}
			summand[j]=term;
	}
	val=KBN_sum(summand,2*n+1);
	val*=(step/3);
	return val;
}

/////////////////////////////////////////////////////////////////////////////////
//////////// \nabla^2 Phi ///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

double fill_arr(double k,double step){
	double val=0;
	double r=0;
	double kr=0;
	for(int j=0;j<(2*n+1)+2;j++){
		r=R_MIN+j*step;
		if(j>1){
			
			val=(sample[j]-2*sample[j-1]+sample[j-2])/(step*step) + (sample[j]-sample[j-2])/(2*step* (r-step));
			
			//val=sample[2*n+1]- sample[j-1];
			kr=k*(r-step);
			
			val*=((r-step)*dbesj0_(&kr));
			
			sample_lap[j-2]=val;	

		}
	}
	
	val=simps_sum(sample_lap,2*n+1,step);
	return(val);
}




double grad_k(double k,double step){
	double val=0;
	double r=0;
	double kr=0;
	for(int j=0;j<(2*n+1)+2;j++){
		r=R_MIN+j*step;
		if(j>1){
			
			val=(sample[j]-2*sample[j-1]+sample[j-2])/(step*step) + (sample[j]-sample[j-2])/(2*step* (r-step));
			
			//val=sample[2*n+1]- sample[j-1];
			kr=k*(r-step);
			
			val*=((r-step)*(r-step)*dbesj1_(&kr));//This line is the only difference
			
			sample_lap[j-2]=val;	

		}
	}
	
	val=simps_sum(sample_lap,2*n+1,step);
	return(val);
}


/////////////////////////////////////////////////////////////////////////////////
//////////// Phi ////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
double fill_arr_2(double k,double step){
	double val=0;
	double r=0;
	double kr=0;
	for(int j=0;j<(2*n+1);j++){
		r=R_MIN+j*step;
		kr=k*(r);
		val=(dbesj0_(&kr)/(r) * sample[j] );
		sample_lap[j]=val;
		//printf("%.3e\n",val);
	}
	
	val=simps_sum(sample_lap,2*n+1,step);
	//printf("%.3e\n",val);
	return(val);
}

double grad_k_2(double k,double step){
	double val=0;
	double r=0;
	double kr=0;
	for(int j=0;j<(2*n+1);j++){
		r=R_MIN+j*step;
		kr=k*(r);
		val=( dbesj0_(&kr)/r - dbesj1_(&kr) *k) * sample[j] ;
		sample_lap[j]=val;	
	}
	
	val=simps_sum(sample_lap,2*n+1,step);
	return(val);
}
