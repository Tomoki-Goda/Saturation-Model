#include<iostream>
#include<cmath>
#include"Kahn.hh"

class Levin{
	private:
		double *a,*s;
		long unsigned int *fact;
		int len;
		int position;
		
		inline double levin_c(int b,int i,int j, int k)const {
			//double val=pow(k+i+1,j-1);
			//val/=pow(k+j+1,j-1);
			double val=((double)(k+i+b))/(k+j+b);
			val=pow(val,j-1);
			return(val);	
		}
		Kahn accum;
	public:
		int reset(){
			position=0;
			Kahn_clear(accum);
			return 0;
		}
		int add_term(double term){
			
			if(position>len){
				printf("Levin:: out of bound\n");
				exit(1);
			}
			a[position]=term;
			accum+=term;
			s[position]=Kahn_total(accum);
			//printf("%.3e %.3e\n",a[position],s[position]);	
			++position;	
			return(position-1);	
		}
		Levin(int len){
			a=(double*)calloc(len+1,sizeof(double));
			s=(double*)calloc(len+1,sizeof(double));
			this->len=len;
			fact=(long unsigned int*)calloc(len,sizeof(long unsigned int));
			fact[0]=1;
			fact[1]=1;
			for(int i=2;i<len;++i){
				fact[i]=i*fact[i-1];
			}
			position=0;
			accum=Kahn_init(3);
		}
		~Levin(){
			free(a);
			free(s);
			Kahn_free(accum);
		}
		inline double sum(int n)const{
			if(n>position-1){
				printf("Levin:: %d st/nd/th element undefined. current=%d\n",n,position-1);
			}
			return(s[n]);
		}
		double accel(int n ,int j)const{
			if(n+j>position-1){
				printf("Levin:: %d st/nd/th element undefined. current=%d\n",n+j+1,position-1);
			}
			double num=0,den=0;
			double val=0;
			Kahn de=Kahn_init(3),nu=Kahn_init(3);
			int b=1;
			for(int i=0;i<=j;++i){
				val=pow(-1,i)*fact[j]/(fact[j-i]*fact[i]);
				val*=levin_c(b,i,j,n) /((b+n+i)*a[n+i]);//Levin type u.
				if(!std::isfinite(val)){
					printf("Levin:: %.3e encountered. a[%d]==%.3e?=0 ??\n",val,n+i,a[n+i]);
					getchar();
				}
				de+=val;	
				nu+=s[n+i]*val;	
			}
			den=Kahn_total(de);
			num=Kahn_total(nu);
			Kahn_free(de);
			Kahn_free(nu);
		    if(den==0.0||!std::isfinite(num)||!std::isfinite(den)){
				printf("Levin:: denominator=%.3e, numerator=%.3e\n",den,num );
			}	
			return(num/den);	
		}
		
};
/*
int main (){
	int n=5;
	int N=15;
	Levin lev(N);
	double levin;
	for(int i=1;i<=N;++i){
		lev.add_term(pow(i,-2));
		
		if(i>=n+2+1){
			levin=lev.accel(i-1-n,n);
			printf("sum = %.3e\t levin = %.6e\t diff = %.3e\t eror=%.3e\n", 
			lev.sum(i-1),levin,lev.sum(i-1)-levin,1.6449340668482264-levin);
		}
	}
}
*/
