#ifndef DCURTIS_H
#define DCURTIS_H
#include<math.h>
#include<stdio.h>
#include<stdlib.h>


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

double dcurtis( double (*func)(const double*), const double min, const double max, const double  eps){
	double smin,smax; //section
	//double x; 
	double total=0, val8=0, val16=0;
	double scale, mid;
	const double prec=1.0e-14;
	double carry=0, carry_total=0;
	double error=0;
	double arg1,arg2;
	const int N=64;//multiple of 4
	double x[N/2-1];
	for(int i=1;i<N/2;i++){
		x[i-1]=cos(i*PI/(N)); 
	}
	smax=max;
	smin=min;
	
	double sample[N/2-1];
	double ends,middle;
	double coeff16[N/2+1],coeff8[N/4+1];
	double val;
	double cos_mat[N/2+1][N/2-1];
	for(int i=0;i<N/2+1;i++){
		for(int n=1;n< N/2;n++){
			cos_mat[i][n-1]=cos(2*n*i*PI/N);
		}
	}

	while(1){
		if((smin>smax)||(smin<min)||(smax>max) ){
			printf("value error [%.3e, %3e] of [%.3e, %.3e]\n",smin,smax,min,max ); 
			getchar();
		}
		scale=(smax-smin)/2;
		mid=(smax+smin)/2;

		if(scale<prec ){
			printf("DCURTIS :: Section size too small. %.3e between %.3e  %.3e of [%.2e, %.2e]\n",scale,smin,smax,min,max );
			getchar();
			return(0);
		}


		for(int i=0;i<N/2-1;i++){
			arg1=mid+x[i]*scale;
			arg2=mid-x[i]*scale;
			sample[i]=(*func)(&arg1)+(*func)(&arg2);
		       	if(isnan(sample[i])==1||isinf(sample[i])==1 ){
				printf("DCURTIS %.5e encountered function(%.3e)+function(%.3e)\n", sample[i],arg1,arg2 );
				printf("%.3e %.3e %.3e %.3e\n",min,max,smin,smax);
				printf("%.3e, %.3e\n",(*func)(&arg1), (*func)(&arg2));
				getchar();
			}	
		}
		
		ends=(*func)(&smin)+(*func)(&smax);
		middle=(*func)(&mid);
		if(isnan(ends*middle)==1||isinf(ends*middle)==1){
			printf("DCURTIS nan  encountered f(1)+f(-1)=%.3e, f(0)=%.3e\n", ends,middle );
			printf("%.3e %.3e %.3e %.3e\n",min,max,smin,smax);
			getchar();
		}	

		val=0;
		for(int i=0;i<N/2+1;i++){
			val=ends/2;
			val+=( ((i/2)*2==i)?(middle):(-middle));
			for(int n=1;n< N/2;n++){
				val+=sample[n-1]*cos_mat[i][n-1];
				if(isnan(val)==1){
					printf("nan sample=%.3e, cos(%.3e) \n", sample[n-1],2*n*i*PI/N);
				}
			}
			coeff16[i]=val;
		}

		val=0;
		for(int i=0;i<N/4+1;i++){
			val=ends/2;
			val+=( ((i/2)*2==i)?(middle):(-middle));
			for(int n=1;n<N/4;n++){
				val+=sample[2*n-1]*cos_mat[i][2*n-1];
			}
			coeff8[i]=val;
		}
		
		val16=coeff16[0]+coeff16[N/2]/(1-N*N);
		for(int i=1;i<N/2;i++){
			val16+=(2*coeff16[i])/(1-4*(i*i));
		}
		val8=coeff8[0]+coeff8[N/4]/(1-N*N/4);
		for(int i=1;i<N/4;i++){
			val8+=(2*coeff8[i])/(1-4*(i*i));
		}


		val16*=2*scale/N;
		val8*=4*scale/N;

		//printf("val16= %.3e  val8= %.3e diff=%.3e [%.3e, %.3e] of [%.3e, %.3e]\n",val16,val8, fabs(val16-val8), smin,smax,min,max);	
		//getchar();
		if(fabs(val16-val8)< eps*(1+fabs(val16) )){
			total=Kahn(total,val16,&carry_total);
			error+=fabs(val16-val8);

			if(fabs(smax-max)<prec){
				//if "end of sector"=="end of integration region", END
				return(total+carry_total);		
			}else{ //move to next sector
				smin=smax;
				smax=max;
			}
		}else{	
			//Half the section
			//smax=mid;
			smax=smin+(scale/2);//quarter

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
		
	return(x*x*exp(-x*x)*sin(PI*x*x));
}

int main(){
	double val=dcurtis(&func , 0,100,1.0e-15);
	printf(" %.5e\n",val);
	return(0);
}
#endif
#endif
