#ifndef DCURTIS_H
#define DCURTIS_H
#include<math.h>
#include<stdio.h>
//#include<stdlib.h>


#define PI 3.1415926535897932384626433
static double Kahn(double a,double b,double *c){
	//double c;
	double	sum=a+b;

	if(a>b){
		*c=sum-a;
		*c=*c-b;
	}else{
		*c=sum-b;
		*c=*c-a;
	}
	return(sum);
}

static const double x16[]={1.000000000000000,	0.980785280403230,	0.923879532511287,	0.831469612302545,	0.707106781186548,	0.555570233019602,	0.382683432365090,	0.195090322016128,	0.000000000000000};	
static const double w16[]={0.062745098039216,	0.298949622697645,	0.603858652345215,	0.871244420655128,	1.111651746945865,	1.305381314253627,	1.451790273891947,	1.540110916903405,	1.571281006575124};	
static const double w8[]={0.126984126984127,	0.584874596864073,	1.117460317460317,	1.446871434881959,	1.574603174603175};	

double dclenshaw(double(*func)(double*),double min,double max,double eps){
	double smin,smax;
	double scale,mid;
	double val16,val8;
	double arg1,arg2;
	double total=0;
	double accum=0;
	const int N=16;
	smin=min;
	smax=max;
	double f[N/2+1];
		
	while(1){
		scale=(smax-smin)/2;
		if(scale<1.0e-15){
			printf("Clenshaw_Curtis::division exceeds limitation. in the domain [%.3e, %.3e] of [%.3e, %.3e] \n",smin,smax,min,max);
			getchar();
		}
		mid=(smax+smin)/2;

		for(int i=0;i<N/2;i++){
			arg1=mid+scale*x16[i];
			arg2=mid-scale*x16[i];
			f[i]=(*func)(&arg1)+(*func)(&arg2);
		}
		f[0]/=2;
		f[N/2]=(*func)(&mid);

		val16=0;
		for(int i=0;i<=N/2;i++){
			val16+=f[i]*w16[i];
			//printf("%.3e\t",f[i]*w[i]);
		}
		//printf("\n");
		val8=0;
		for(int i=0;i<=N/4;i++){
			val8+=f[2*i]*w8[i];
		}
		val16*=2*scale/N;
		val8*=4*scale/N;
		//printf("%.5e , %.5e in the domain [%.3e, %.3e] of [%.3e, %.3e] \n",valfull,valhalf,smin,smax,min,max);
		if(fabs(val16-val8)<eps*(1+fabs(val16)) ){
			//total+=valfull;
			total=Kahn(total,val16,&accum);
			
			if(fabs(smax-max)<1.0e-15){
			return(total+accum);
			//bereak;
			}
			smin=smax;
			smax=max;
		}else{
			//smax=mid;
			smax=smin+(scale/2);
		}
	}
}

//////////////////////////////////////////////////////////
//
//
//
#define TEST 0
#if TEST==1
double func(double *X){
	double x=*X;
		
	return(x*x*exp(-x*x));
}

int main(){
	double val=dclenshaw(&func , 0,100,1.0e-15);
	printf(" %.5e\n",val);
	return(0);
}
#endif
#endif
