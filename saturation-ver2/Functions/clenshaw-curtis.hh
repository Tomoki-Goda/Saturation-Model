#include<cmath>
#include<iostream>
//#define PI 3.141592653589793238462643383279502884197
#include<quadmath.h>


//PI 3.141592653589793238462643383279502884197


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
		//__float128 c;
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
		__float128 quad_Kahn(__float128 a,__float128 b,__float128 *c){
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

		int prepare(int N,__float128* x,__float128* w){

			x[0]=1;
			x[1]=cosq(PI/N);
			for(int i =2;i<N/2;i++){
				//x[i]=2*x[1]*x[i-1]-x[i-2];
				x[i]=cosq(i*PI/N);
			}
			x[N/2]=0;
			
			char str[200];
			for(int i=0;i<N/2+1;i++){
				//str=(char*)malloc(100);
				quadmath_snprintf(str,sizeof(str),"%.30Qe",x[i]);
				fprintf(file,"%s, ",str);	
				//printf("%s",str);
				//printf("%.20e\n", (double)x[i]);	
				//free(str);
			}
			fprintf(file,"\n");
			__float128 c[N/2+1];
			c[0]=1;
			c[N/2]=1.0/(1-N*N);
			for(int i=1;i<N/2;i++){
				c[i]=2.0/(1-4*i*i);
			}
			//for(int i=0;i<N/2+1;i++){	
			//	printf("%.3e\t",c[i]);
			//}
			//printf("\n\n");
			__float128 t[N/2+1];
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
			}
			
			for(int i=0;i<N/2+1;i++){
				//str=(char*)malloc(100);
				quadmath_snprintf(str,sizeof(str),"%.30Qe",w[i]);
				fprintf(file,"%s, ",str);	
				//printf("%s",str);	
				//fprintf(file,"%.25Qe,\t",w[i]);
				//free(str);
			}
			//free(str);
			fprintf(file,"\n");
			return 0;
		}
		int N=128;
		__float128 x[128/2+1];
		__float128 wfull[128/2+1];
		__float128 xhalf[128/4+1];
		__float128 whalf[128/4+1];


		__float128 f[128/2+1];

	public:
		int Print=0;
		int DIV=2;
		double  integrate(double(*func)(double*),double min,double max,double eps){
			double smin,smax;
			double scale,mid;
			double valfull,valhalf;
			double arg1,arg2;
			double total=0;
			double accum=0,accum2=0;

			smin=min;
			smax=max;
			int licz=0;
			int totlicz=0;
			double efficiency;
			while(true){
				totlicz++;

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
				accum2=0;
				for(int i=0;i<=N/2;i++){
					valfull=Kahn(valfull,f[i]*wfull[i],&accum2);
					//printf("%.3e\t",f[i]*w[i]);
				}
				valfull+=accum2;
				//printf("\n");
				valhalf=0;
				accum2=0;
				for(int i=0;i<=N/4;i++){
					valhalf=Kahn(valhalf,f[2*i]*whalf[i],&accum2);
				}
				valhalf+=accum2;

				valfull*=2*scale/N;
				valhalf*=4*scale/N;
				//printf("%.5e , %.5e in the domain [%.3e, %.3e] of [%.3e, %.3e] \n",valfull,valhalf,smin,smax,min,max);
				if(fabs(valfull-valhalf)<eps*(1+fabs(valfull)) ){
					//total+=valfull;
					total=Kahn(total,valfull,&accum);
					licz++;
					if(fabs(smax-max)<1.0e-15){
						efficiency=((double)licz)/totlicz;
						if(Print==1||efficiency<0.5){
							//DIV*=2;
							printf("Clenshaw_Curtis:: Efficiency %d/%d=%f\n",N*licz,N*totlicz,((double)licz)/totlicz);
						}
						return(total+accum);
						//bereak;
					}
					smin=smax;
					//smax=max;
					smax=( ( (max-(2*DIV*scale+smin))<1.0e-10)?(max):(2*DIV*scale+smin) );
				}else{
					//smax=mid;
					smax=smin+(scale/DIV);
				}
			}
		}

};
/*
double func(double *X){
	double x=*X;
	double val=(1-50*x*x)*exp(-x*x);
	return(val);
}
int main(){
	Clenshaw_Curtis integration(16);
	double res=integration.integrate(&func,0,100,1.0e-10);
	std::cout<<res<<std::endl;
 	return(0);
	
}*/
