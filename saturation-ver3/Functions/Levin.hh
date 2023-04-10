#include<iostream>
#include<cmath>
#include"Kahn.hh"

class Levin{
	private:
		double *a,*s;
		long unsigned int *fact;
		//long unsigned int len;
		long unsigned int position;
		
		inline double levin_c(int b,int i,int j, int k)const {
			double val=((double)(k+i+b))/(k+j+b);
			val=pow(val,j-1);
			return(val);	
		}
		Kahn accum;
		int clen=0;//length of record
	public:
		int reset(){
			position=0;
			Kahn_clear(accum);
			return 0;
		}
		int add_term(double term){
			a[position-clen*(position/clen)]=term;
			accum+=term;
			s[position-clen*(position/clen)]=Kahn_total(accum);
			++position;	
			return(position-1);	
		}
		Levin(int len){
			clen=len;
			fact=(long unsigned int*)calloc(clen,sizeof(long unsigned int));
			a=(double*)calloc(clen,sizeof(double));
			s=(double*)calloc(clen,sizeof(double));
			//this->len=len;
			fact[0]=1;
			fact[1]=1;
			for(int i=2;i<clen;++i){
				fact[i]=i*fact[i-1];
			}
			position=0;
			accum=Kahn_init(3);
		}
		~Levin(){
			free(a);
			free(s);
			free(fact);
			Kahn_free(accum);
		}
		inline double sum(int n)const{
			
			if(n>position-1||n+clen<position){
				printf("Levin::sum %d st/nd/th element undefined. current=%lu\n",n,position-1);
			}
			return(s[n-clen*(n/clen)]);
		}
		double accel(int n ,int j)const{
			
			if(n+j>position-1||n+clen<position){
				printf("Levin::accel %d st/nd/th (%d elements from %d to %d of %d) element undefined. current=%lu\n",n+j,j,n+j-clen*((n+j)/clen),n-clen*((n)/clen),clen, position-1);
				printf("%d %lu %d %lu \n",n+j,position-1,n,position-clen);
				return -999;
			}
			double num=0,den=0;
			double val=0;
			Kahn de=Kahn_init(3),nu=Kahn_init(3);
			int b=1,cyc=0;
			
			for(int i=0;i<=j;++i){
				cyc=n+i-clen*((n+i)/clen);
				//n+i th position corresponds to cyc th position in the record
				
				val=pow(-1,i)*fact[j]/(fact[j-i]*fact[i]);
				val*=levin_c(b,i,j,n) /((b+n+i)*a[cyc]);//Levin type u.
				
				if(!std::isfinite(val)){
					printf("Levin:: %.3e encountered.\n",val);
					for(int j=0;j<clen;++j){
						printf(" a[%d]==%.3e \t s[%d]==%.3e ??\n",n+i-j,a[n+i-j-clen*((n+i-j)/clen)],n+i-j,s[n+i-j-clen*((n+i-j)/clen)]);
					}
					printf("\n");
					getchar();
				}
				de+=val;	
				nu+=s[cyc]*val;	
				
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
