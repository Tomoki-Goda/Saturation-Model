#include<cmath>
#include<iostream>
#ifndef PI
	#define PI 3.1415926535897932384626433
#endif
class Clenshaw_Curtis{
	public:
	Clenshaw_Curtis(int n){
		file=fopen("x-and-w.txt","w");	
		N=4*(n/4);
		if(N!=n){
			printf("using N=%d\n",N);
		}
		if(N>128){
			printf("max 128\n");
			N=128;
		}
		prepare(N,x,wfull);
		prepare(N/2,xhalf,whalf);
		fclose(file);
	}
	private:
		FILE * file;
		double Kahn(double a,double b,double *c){
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


		int prepare(int N,double* x,double* w){

			x[0]=1;
			x[1]=cos(PI/N);
			for(int i =2;i<(N/2+1);i++){
				//x[i]=2*x[1]*x[i-1]-x[i-2];
				x[i]=cos(i*PI/N);
			}
			double c[N/2+1];
			c[0]=1;
			c[N/2]=1.0/(1-N*N);
			for(int i=1;i<N/2;i++){
				c[i]=2.0/(1-4*i*i);
			}
			double t[N/2+1];
			t[0]=1;
			for(int i=0;i<N/2+1;i++){
				t[1]=cos(2*i*PI/N);
				w[i]=c[0]+c[1]*t[1];
				for(int j=2;j<N/2+1;j++){
					//t[j]=2*t[1]*t[j-1]-t[j-2];
					t[j]=cos(i*j*2*PI/N);
					w[i]+=t[j]*c[j];
				}
			}
			return 0;
		}
		int N=128;
		double x[128/2+1];
		double wfull[128/2+1];
		double xhalf[128/4+1];
		double whalf[128/4+1];


		double f[128/2+1];
	public:
		double  integrate(double(*func)(double*),double min,double max,double eps){
			double smin,smax;
			double scale,mid;
			double valfull,valhalf;
			double arg1,arg2;
			double total=0;
			double accum=0;

			smin=min;
			smax=max;
			
			while(true){
				scale=(smax-smin)/2;
				if(scale<1.0e-15){
					printf("Clenshaw_Curtis::division exceeds limitation. in the domain [%.3e, %.3e] of [%.3e, %.3e] \n",smin,smax,min,max);
					getchar();
				}
				mid=(smax+smin)/2;

				for(int i=0;i<N/2;i++){
					arg1=mid+scale*x[i];
					arg2=mid-scale*x[i];
					f[i]=(*func)(&arg1)+(*func)(&arg2);
				}
				f[0]/=2;
				f[N/2]=(*func)(&mid);

				valfull=0;
				for(int i=0;i<=N/2;i++){
					valfull+=f[i]*wfull[i];
					//printf("%.3e\t",f[i]*w[i]);
				}
				//printf("\n");
				valhalf=0;
				for(int i=0;i<=N/4;i++){
					valhalf+=f[2*i]*whalf[i];
				}
				valfull*=2*scale/N;
				valhalf*=4*scale/N;
				//printf("%.5e , %.5e in the domain [%.3e, %.3e] of [%.3e, %.3e] \n",valfull,valhalf,smin,smax,min,max);
				if(fabs(valfull-valhalf)<eps*(1+fabs(valfull)) ){
					//total+=valfull;
					total=Kahn(total,valfull,&accum);
					
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

};
/*
double func(double *X){
	double x=*X;
	double val=x*x*exp(-x*x);
	return(val);
}
int main(){
	Clenshaw_Curtis integration(16);
	double res=integration.integrate(&func,0,100,1.0e-10);
	std::cout<<res<<std::endl;
 	return(0);
	
}
*/
