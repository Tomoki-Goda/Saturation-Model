#ifndef DCLENSHAW_H
#define DCLENSHAW_H
#include<math.h>
#include<stdio.h>
#include"./Kahn.h"
//#include<stdlib.h>

#ifndef PI
#define PI 3.141592653589793238462643383279502884197
#endif


static const double x16[]={1.000000000000000000000000000000e+00, 9.807852804032304506194118511447e-01, 9.238795325112867619863336960547e-01, 8.314696123025452498358628945367e-01, 7.071067811865475460497457679922e-01, 5.555702330196022565633495419049e-01, 3.826834323650898141569391948192e-01, 1.950903220161283203970903689557e-01, 0.000000000000000000000000000000e+00};	
static const double w16[]={6.274509803921572703711007079619e-02, 2.989496226976449054167638106882e-01, 6.038586523452146797285979221456e-01, 8.712444206551274286632311850505e-01, 1.111651746945864560764091023555e+00, 1.305381314253626384086121496437e+00, 1.451790273891946807960151513108e+00, 1.540110916903405124741484817321e+00, 1.571281006575124188778813660861e+00 };	
static const double w8[]={1.269841269841270256502063773496e-01, 5.848745968640726177649829569567e-01, 1.117460317460317429023630586917e+00, 1.446871434881959070257763581104e+00, 1.574603174603174567114383108901e+00};


template<typename TYPE,typename args_type>static double dclenshaw(const TYPE &func,const args_type par , const double a, const double b, const double eps){
	int MAX_RECURSION=15;
//double dclenshaw(double(*func)(double*),double a,double b,double eps){
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
	double valfull,valhalf;
	double arg1,arg2;
	double total=0;
	double accum[3]={0};
	double accum2[3]={0};
	const int N=16;
	double increase;
	smin=min;
	smax=max;

	int licz=0,licztot=0 , counter=0;
	double f[N/2+1];
	if(fabs(min-max)<1.0e-15){
		return(0);
	}		
	while(1){
		counter++;
		licztot++;
		scale=(smax-smin)/2;
		mid=(smax+smin)/2;

		for(int i=0;i<N/2;i++){
			arg1=mid+scale*x16[i];
			arg2=mid-scale*x16[i];
			f[i]=func(&arg1,par)+func(&arg2,par);
			if((isnan(f[i])+isinf(f[i]))!=0){
				printf("%.3e encountered at %.3e or %.3e\n",f[i],arg1,arg2 );
			}
		}
		f[0]/=2;
		f[N/2]=func(&mid,par);

		if((isnan(f[N/2])+isinf(f[N/2]))!=0){
			printf("%.3e encountered at %.3e \n",f[N/2],mid );
		}
		valfull=0;
		Kahn_init(accum2,3);
		for(int i=0;i<=N/2;i++){
			valfull=Kahn_Sum(valfull,f[i]*w16[i],accum2,3);
		}
		valfull=Kahn_total(valfull,accum2,3);
		valhalf=0;
		Kahn_init(accum2,3);
		for(int i=0;i<=N/4;i++){
			valhalf=Kahn_Sum(valhalf,f[2*i]*w8[i],accum2,3);
		}
		//valhalf+=accum2;
		valhalf=Kahn_total(valhalf,accum2,3);

		valfull*=2*scale/N;
		valhalf*=4*scale/N;
		if(( fabs(valfull-valhalf)<eps*(fabs(valfull)) ) || (  fabs(valfull-valhalf)<eps )|| (counter==MAX_RECURSION)){//Need improvement
			if(counter==MAX_RECURSION){
				printf("MAX_RECURSION:: evaluated %d times. Increase MAX_RECURSION\n",MAX_RECURSION );
				printf("[%.3e, %.3e] of [%.3e, %.3e] after %d / %d \n",smin,smax,min,max, licz,licztot);
				printf("valfull= %.3e , valhalf= %.3e  diff=%.3e\n",valfull,valhalf,valfull-valhalf);
				for(int i=0;i<N/2;i++){
					arg2=mid-scale*x16[i];
					printf("%.3e\n",func(&arg2,par));
				}
				arg2=mid;
				printf("%.3e\n",func(&arg2,par));
				for(int i=0;i<N/2;i++){
					arg2=mid+scale*x16[N/2-i-1];
					printf("%.3e\n",func(&arg2,par));
				}
				printf("\n");
			}
			
			total=Kahn_Sum(total,valfull,accum,3);
			licz++;
			counter=0;
			if(fabs(smax-max)<1.0e-15){
				return(sign*Kahn_total(total,accum,3) );
			}
			smin=smax;
			increase=(4*scale);
			smax=((max-(smin+increase)<(increase/2))?(max):(smin+increase));
		}else{
			smax=smin+(scale/2);
		}
		if((smax-smin)<1.0e-15){
			printf("Abs: %.5e %.5e \tRel: %.5e %.5e\n ",fabs(2*valfull*(max-min)/scale),eps,fabs(valfull-valhalf),eps*(fabs(valfull)));
			printf("Clenshaw_Curtis::division exceeds limitation. in the domain [%.3e, %.3e] of [%.3e, %.3e] after %d / %d \n",smin,smax,min,max, licz,licztot);						
			printf("valfull= %.3e , valhalf= %.3e \n",valfull,valhalf);
			//getchar();
		}
	}
}


#endif
