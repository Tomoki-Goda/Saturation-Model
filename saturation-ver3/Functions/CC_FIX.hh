#ifndef CC_FIX_LOADED
#define CC_FIX_LOADED 1

#include<cmath>
#include<iostream>
#ifndef PI
#define PI 3.141592653589793238462643383279502884197
#endif
#include<quadmath.h>
#include"./kahnsum.hh"
//#include"./Kahn.h"
//#include"./clenshaw.h"


//PI 3.141592653589793238462643383279502884197

class Clenshaw_Curtis_FIX{
	
	public:
		unsigned int Print=0;
		unsigned int DIV=2;
	protected:
		
		int N=128;
		__float128 *__restrict__ x;//=(__float128*)malloc((N/2+1)*sizeof(__float128) );
		__float128 *__restrict__ wfull;
		__float128 *__restrict__ xhalf;
		__float128 *__restrict__ whalf;
		//__float128 x[128/2+1];
		//__float128 wfull[128/2+1];
		//__float128 xhalf[128/4+1];
		//__float128 whalf[128/4+1];
	//int N;
	//__float128 *x;
	//__float128 *wfull;
	//__float128 *xhalf;
	//__float128 *whalf;
	public:
	
	
	Clenshaw_Curtis_FIX(int n){
		//file=fopen("x-and-w.txt","w");	
		N=4*(n/4);
		if(N!=n){
			printf("using N=%d\n",N);
		}
		if(N>128){
			printf("max 128\n");
			N=128;
		}
		
		//
		x=(__float128*)malloc((n/2+1)*sizeof(__float128) );
		wfull=(__float128*)malloc((n/2+1)*sizeof(__float128) );
		xhalf=(__float128*)malloc((n/4+1)*sizeof(__float128) );
		whalf=(__float128*)malloc((n/4+1)*sizeof(__float128) );
		prepare(N,x,wfull);
		prepare(N/2,xhalf,whalf);
	}
	~Clenshaw_Curtis_FIX(){
		
	free(x);
	free(wfull);
	free(xhalf);
	free(whalf);	
	}
	private:
		__float128 quad_Kahn(__float128 a,__float128 b,__float128 *c)const{
		//__float128 c;
			__float128	sum=a+b;

			if(a>b){
				*c=sum-a;
				*c=*c-b;
			}else{
				*c=sum-b;
				*c=*c-a;
			}
			return(sum);
		}

		int prepare(int N,__float128* __restrict__ x,__float128* __restrict__ w)const{

			x[0]=1;
			x[1]=cosq(PI/N);
			for(int i =2;i<N/2;i++){
				x[i]=cosq(i*PI/N);
			}
			x[N/2]=0;
			
			//char str[200];
			//for(int i=0;i<N/2+1;i++){
			//	quadmath_snprintf(str,sizeof(str),"%.30Qe",x[i]);
			//}
			//fprintf(file,"\n");
			__float128 *__restrict__ c=(__float128*)malloc((N/2+1)*sizeof(__float128));
			c[0]=1;
			c[N/2]=1.0/(1-N*N);
			for(int i=1;i<N/2;i++){
				c[i]=2.0/(1-4*i*i);
			}
			//__float128 t[N/2+1];
			__float128 *__restrict__ t=(__float128*)malloc((N/2+1)*sizeof(__float128));
			__float128 accum=0;
			t[0]=1;
			for(int i=0;i<N/2+1;i++){
				t[1]=cosq(2*i*PI/N);
				w[i]=c[0]+c[1]*t[1];
				for(int j=2;j<N/2+1;j++){
					//t[j]=2*t[1]*t[j-1]-t[j-2];
					t[j]=cosq(i*j*2*PI/N);
					w[i]=quad_Kahn(w[i],t[j]*c[j],&accum);
				}
				w[i]+=accum;
				if(i==0){
				w[i]/=2;
				}
			}
			
			//for(int i=0;i<N/2+1;i++){
			//	quadmath_snprintf(str,sizeof(str),"%.30Qe",w[i]);
			//}
			//free(str);
			free(c);
			free(t);
			//fprintf(file,"\n");
			return 0;
		}
		

	protected:
				
		double integral16(
				double(*func)(double*,void*),
				void* __restrict__ args,double smin,double smax,
			       	double * __restrict__ res, double* __restrict__ half )const{
			double* __restrict__ f=(double*)malloc((N/2+1)*sizeof(double));
			//double ( __restrict__ f)[N/2+1];
			double arg1=0,arg2=0,scale=0,mid=0;
			double valfull=0,valhalf=0;
			Kahn total_sect(3);			
			total_sect.zero();
			
			scale=(smax-smin)/2;

			mid=(smax+smin)/2;
			for(int i=0;i<N/2;i++){
				arg1=mid+scale*x[i];
				arg2=mid-scale*x[i];
				f[i]=(*func)(&arg1,args)+(*func)(&arg2,args);
				if((std::isnan(f[i])+std::isinf(f[i]))!=0){
					printf("%f encountered at %e or %e of [%.3e, %.3e]\n",f[i],arg1,arg2, smin,smax);
					//getchar();
				}
			}
			f[N/2]=(*func)(&mid,args);
			if((std::isnan((f[0]*f[N/2]))+std::isinf((f[0]*f[N/2])))!=0){
				printf("%f encountered at %e of [%.3e, %.3e]\n", f[N/2],mid,smin,smax);
				//getchar();
			}
			
			for(int i=0;i<=N/2;i++){
				total_sect+=f[i]*wfull[i];
				//total_sect=Kahn(total_sect,f[i]*wfull[i],accum,3);
			}
			//valfull=Kahn_total(total_sect,accum,3);
			valfull=total_sect.total(0);

			for(int i=0;i<=N/4;i++){
				total_sect+=f[2*i]*whalf[i];
				//total_sect=Kahn(total_sect,f[2*i]*whalf[i],accum,3);
			}
			valhalf=total_sect.total(0);
			valfull*=2*scale/N;
			valhalf*=4*scale/N;
			*res=valfull;
			*half=valhalf;
			free(f);
			return(valfull);
		}
};
#endif
