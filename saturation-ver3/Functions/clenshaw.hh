#ifndef DCLENSHAW_HH
#define DCLENSHAW_HH 1
#include<math.h>
#include<stdio.h>
#include<string>
#include"./Kahn.h"
//#include<stdlib.h>

#ifndef PI
#define PI 3.141592653589793238462643383279502884197
#endif
typedef  struct{int N=128; double wfull[65], whalf[33], x[65]; std::string tag="unnamed"; int max_rec=7;} CCIntegral;
template<typename TYPE,typename args_type>static double dclenshaw(const CCIntegral &data,TYPE &func, args_type par , const double a, const double b, const double eps, const double Aeps){
	const double (&x16)[]=data.x;
	const double (&w16)[]=data.wfull;
	const double (&w8)[]=data.whalf;
	const int N=data.N;
	int MAX_RECURSION=data.max_rec;
//double dclenshaw(double(&func)(const double*, const void*),const void* args,const double a,const double b,const double eps){
	double sign, max,min;

	if(b>a){
		min=a;
		max=b;
		sign=1;
	}else if(a>b){
		max=a;
		min=b;
		sign=-1;
	}else{
		return(0);
	}
	double smin,smax;
	double scale,mid;
	double valfull,valhalf;
	double arg;
	double total=0;
	double accum[3]={0};
	double accum2[3]={0};
	//const int N=16;
	double increase;
	smin=min;
	//smax=max/data.div;
	smax=max;//data.init_div;
	int licz=0,licztot=0 , counter=0;
	double f[N/2+1];
	if(fabs(min-max)<0.0){
		return(0);
	}		
	while(1){
		++counter;
		++licztot;
		scale=(smax-smin)/2;
		mid=(smax+smin)/2;

		for(int i=0;i<N/2;i++){
			arg=scale*x16[i];
			f[i]=func(mid+arg,par)+func(mid-arg,par);
		}
		
		f[N/2]=func(mid,par);
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
		valhalf=Kahn_total(valhalf,accum2,3);

		valfull*=2*scale/N;
		valhalf*=4*scale/N;
#if DCLENSHAW_HH==1		
		if(not(std::isfinite(valfull)&&std::isfinite(valhalf))){
			printf("Clenshaw_Curtis:: in \"%s\" %.3e  %.3e encountered\n",(data.tag).c_str(),valfull,valhalf);
			goto Error;
		}
#endif
		
		
		if(( fabs(valfull-valhalf)<eps*(fabs(valfull)) ) || (  fabs(valfull-valhalf)<Aeps )|| (counter==MAX_RECURSION)){//Need improvement
			if(counter==MAX_RECURSION){
				printf("Clenshaw_Curtis::MAX_RECURSION:: in \"%s\", evaluated %d times. Increase MAX_RECURSION\n",(data.tag).c_str(), MAX_RECURSION );
				printf("[%.3e, %.3e] of [%.3e, %.3e] after %d / %d \n",smin,smax,min,max, licz,licztot);
				printf("valfull= %.3e , valhalf= %.3e  diff=%.3e\n",valfull,valhalf,valfull-valhalf);
				//getchar();
				goto Error;
			}
			
			total=Kahn_Sum(total,valfull,accum,3);
			++licz;
			//data.div+=pow(4,counter);
			
			counter=0;
			if(fabs(smax-max)==0.0){
				//data.div/=licz;
				return(sign*Kahn_total(total,accum,3) );
			}
			smin=smax;
			increase=(4*scale);
			smax=((max-(smin+increase)<(increase/2))?(max):(smin+increase));
		}else{
			smax=smin+(scale/2);
		}
		if(((max-min)-(smax-smin))==(max-min)){
			printf("Clenshaw_Curtis:: in \"%s\", division exceeds limitation. in the domain [%.3e, %.3e] of [%.3e, %.3e] after %d / %d \n",(data.tag).c_str(), smin,smax,min,max, licz,licztot);						
			printf("valfull= %.3e , valhalf= %.3e \n",valfull,valhalf);
			//getchar();
			goto Error;
		}
		
	}
	
	Error:
		for(int i=0;i<N/2;i++){
			arg=mid-scale*x16[i];
			printf("f(%.3e) = %.3e\n",arg,func(arg,par));
		}
		arg=mid;
		printf("f(%.3e) = %.3e\n",arg,func(arg,par));
		for(int i=0;i<N/2;i++){
			arg=mid+scale*x16[N/2-i-1];
			printf("f(%.3e) = %.3e\n",arg, func(arg,par));
		}
		printf("\n");
		getchar();
		return 0;
	
}



CCIntegral CCprepare(const int N){
	CCIntegral data;
	
	if((N/4)*4!=N || N>128){
		printf("N=%d has to be multiple of 4, <= 128\n",N );
	}
	data.N=N;
	data.x[0]=1;
	data.x[1]=cos(PI/N);
	for(int i =2;i<N/2;i++){
		data.x[i]=cos(i*PI/N);
	}
	data.x[N/2]=0;
	
	double *__restrict__ c=(double*)malloc((N/2+1)*sizeof(double));
	c[0]=1;
	c[N/2]=1.0/(1-N*N);
	for(int i=1;i<N/2;i++){
		c[i]=2.0/(1-4*i*i);
	}
	double *__restrict__ t=(double*)malloc((N/2+1)*(N/2+1)*sizeof(double));
	double accum[3]={0};
	for(int i=0;i<N/2+1;i++){
		for(int j=0;j<N/2+1;j++){
			if(( (i*j)/N)*N==i*j){
				t[i*(N/2+1)+j]=1;
			}else if(( (2*i*j)/N)*N==2*i*j){
				t[i*(N/2+1)+j]=-1;
			}else if(( (4*i*j)/N)*N==4*i*j){
				t[i*(N/2+1)+j]=0;
			}else{
				t[i*(N/2+1)+j]=cos(2*i*j*PI/N);
			}
			//t[i*(N/2+1)+j]=((i*j==0)?(1):(cos(2*i*j*PI/N)) );
		}
	}
	
	for(int i=0;i<N/2+1;i++){
		Kahn_init(accum,3);
		//t[1]=cos(2*i*PI/N);
		data.wfull[i]=c[0]+c[1]*t[i*(N/2+1)+1];
		
		for(int j=2;j<N/2+1;j++){
			//t[j]=cos(i*j*2*PI/N);
			data.wfull[i]=Kahn_Sum(data.wfull[i],t[i*(N/2+1)+j]*c[j],accum,3);
		}
		data.wfull[i]=Kahn_total(data.wfull[i],accum,3);
		if(i==0){
			data.wfull[i]/=2;
		}
		//printf("%.10e\t",data.wfull[i]);
	}//printf("\n");
	
	//t[0]=1;
	for(int i=0;i<N/4+1;i++){
		Kahn_init(accum,3);
		//t[1]=cos(4*i*PI/N);
		data.whalf[i]=c[0]+c[1]*t[i*(N/2+1)+2];
		for(int j=2;j<N/4+1;j++){
			//t[j]=cos(i*j*4*PI/N);
			data.whalf[i]=Kahn_Sum(data.whalf[i],t[i*(N/2+1)+2*j]*c[j],accum,3);
		}
		data.whalf[i]=Kahn_total(data.whalf[i],accum,3);
		
		//data.whalf[i]+=accum;
		if(i==0){
			data.whalf[i]/=2;
		}
		//printf("%.10e\t",data.whalf[i]);
	}//printf("\n");
	
	free(c);
	free(t);
	//getchar();
	return data;
}
CCIntegral CCprepare(const int N,const std::string &tag){
	CCIntegral data=CCprepare(N);
	data.tag=tag;
	
	return data;
}
CCIntegral CCprepare(const int N,const std::string &tag,int d){
	CCIntegral data=CCprepare(N);
	data.tag=tag;
	data.max_rec=d;
	
	return data;
}


#endif
