#ifndef DCURTIS_H
#define DCURTIS_H
#include<math.h>
#include<stdio.h>
//#include<stdlib.h>


#define PI 3.141592653589793238462643383279502884197
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

static const double x16[]={1.000000000000000000000000000000e+00, 9.807852804032304506194118511447e-01, 9.238795325112867619863336960547e-01, 8.314696123025452498358628945367e-01, 7.071067811865475460497457679922e-01, 5.555702330196022565633495419049e-01, 3.826834323650898141569391948192e-01, 1.950903220161283203970903689557e-01, 0.000000000000000000000000000000e+00};	
static const double w16[]={6.274509803921572703711007079619e-02, 2.989496226976449054167638106882e-01, 6.038586523452146797285979221456e-01, 8.712444206551274286632311850505e-01, 1.111651746945864560764091023555e+00, 1.305381314253626384086121496437e+00, 1.451790273891946807960151513108e+00, 1.540110916903405124741484817321e+00, 1.571281006575124188778813660861e+00 };	
static const double w8[]={1.269841269841270256502063773496e-01, 5.848745968640726177649829569567e-01, 1.117460317460317429023630586917e+00, 1.446871434881959070257763581104e+00, 1.574603174603174567114383108901e+00};
double dclenshaw(double(*func)(double*),double a,double b,double eps){
	double sign, max,min;

	if((1-(a-b))==1){
		return(0);
	}else if(a>b){
		max=a;
		min=b;
		sign=-1;
	}else{
		min=a;
		max=b;
		sign=1;
	}
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
	if(fabs(min-max)<1.0e-15){
		return(0);
	}		
	while(1){
		scale=(smax-smin)/2;
		if(scale<2.0e-15){
			printf("Clenshaw_Curtis::division exceeds limitation. in the domain [%.3e, %.3e] of [%.3e, %.3e] scale = %.5e\n",smin,smax,min,max,scale);
			getchar();
		}
		mid=(smax+smin)/2;

		for(int i=0;i<N/2;i++){
			arg1=mid+scale*x16[i];
			arg2=mid-scale*x16[i];
			f[i]=(*func)(&arg1)+(*func)(&arg2);
			if((isnan(f[i])+isinf(f[i]))!=0){
				printf("%.3e encountered at %.3e or %.3e\n",f[i],arg1,arg2 );
			}
		}
		f[0]/=2;
		f[N/2]=(*func)(&mid);

		if((isnan(f[0]*f[N/2])+isinf(f[0]*f[N/2]))!=0){
			printf("%.3e encountered at %.3e \n",f[N/2],mid );
		}
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
				return(sign*(total+accum) );
			//bereak;
			}
			smin=smax;
			//smax=max;
			smax=((max-(smin+2*4*scale)<1.0e-10)?(smin+2*4*scale):(max));
		}else{
			//smax=mid;
			smax=smin+(scale/4);
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
