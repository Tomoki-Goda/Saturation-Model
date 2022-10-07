
//extern double k_group_sum(const double *arr,int len);
//#include"./kahnsum.h"
//#include"./critical-line.h"

#include"../Functions/kahnsum.h"
//#include"./main.h"
extern double dbesj0_(double*);
extern double dbesj1_(double*);
extern double mod_x(double, double,int);
extern double SIGMA(double , double ,double , const double *, const double*);
extern double BASE_SIGMA(double , double ,double , const double *);

extern int parameter(const double*,double*,double*);
extern void approx_xg(const double *);



static const int n=500;
static double sample[2*500+1 +2]={0};
static double sample_lap[2*500+1]={0};
static double summand[2*500+1]={0};

double sample_sigma(double * sample , double step, double x,double Q2,const double * sigpar, const double * sudpar){ 
	// one only needs to compute a set of sigma once.
	double r=0;
	double xm=0;
	double val=0;
	for(int j=0;j<(2*n+1)+2;j++){
		r=R_MIN+j*step;
		
		//val=0;
		//for(int i=0; i<(NF-1); i++){
		//for(int i=0; i<1; i++){
		//xm=mod_x(x,Q2,0);
#if ((MODEL==0||MODEL==1)||(PHI==1))
		val=SIGMA(r,x,Q2,sigpar,sudpar);
#else 
		val=BASE_SIGMA(r,x,Q2,sigpar);
#endif
			
			//printf("val=%.3e %.3e %.3e %.3e\n",val,r ,xm,Q2);
		//}
		//printf("%.3e\n",val);
		sample[j]=val;
	}
}

double simps_sum(double * sample, int len ,  double step){
	double summand[len];
	double term, val;
	for(int j=0;j<(2*n+1);j++){	
			term=sample[j];
				
			if((j==0)||(j==2*n)){
			
			} else if( (j/2)*2==j ){
				term*=2;
			}
			else{
				term*=4;	
			}
			summand[j]=term;
	}
	val=k_group_sum(summand,2*n+1);
	val*=(step/3);
	return val;
}

/////////////////////////////////////////////////////////////////////////////////
//////////// \nabla^2 Phi ///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
extern double exp_sud(double,double , double);
extern int compute_mu2(double,const double *,const double *,int);

double fill_arr(double k,double step,double *sudpar,double q2){
	double val=0;
	double r=0;
	double kr=0;
	double mu2_arr[1]={0};
	int signal;
	for(int j=0;j<(2*n+1)+2;j++){
		r=R_MIN+j*step;
		if(j>1){
			signal=compute_mu2(r, sudpar, mu2_arr,1);//compute mu2
			
			val=(sample[j]-2*sample[j-1]+sample[j-2])/(step*step) + (sample[j]-sample[j-2])/(2*step* (r-step));
			
			val*=exp_sud(r,mu2_arr[0],q2);
			//val=sample[2*n+1]- sample[j-1];
			kr=k*(r-step);
			val*=((r-step)*dbesj0_(&kr));
			sample_lap[j-2]=val;
		}
	}
	
	val=simps_sum(sample_lap,2*n+1,step);
	//printf("%.5e %.5e %.5e\n",q2, k*k, val);
	return(val);
}




double grad_k(double k,double step,double *sudpar,double q2){
	double val=0;
	double r=0;
	double kr=0;
	double mu2_arr[1]={0};
	int signal;
	for(int j=0;j<(2*n+1)+2;j++){
		r=R_MIN+j*step;
		if(j>1){
			signal=compute_mu2(r, sudpar, mu2_arr,1);//compute mu2

			val=(sample[j]-2*sample[j-1]+sample[j-2])/(step*step) + (sample[j]-sample[j-2])/(2*step* (r-step));
			val*=exp_sud(r,mu2_arr[0],q2);
			//val=sample[2*n+1]- sample[j-1];
			kr=k*(r-step);
			
			val*=((r-step)*(r-step)*dbesj1_(&kr));//This line is the only difference
			
			sample_lap[j-2]=val;	

		}
	}
	
	val=simps_sum(sample_lap,2*n+1,step);
	return(val);
}



/////////////////////////////////////////////////////////////////////////////////
//////////// Phi ////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
double fill_arr_2(double k,double step){
	double val=0;
	double r=0;
	double kr=0;
	for(int j=0;j<(2*n+1);j++){
		r=R_MIN+j*step;
		kr=k*(r);
		val=(dbesj0_(&kr)/(r) * sample[j] );
		sample_lap[j]=val;
		//printf("%.3e\n",val);
	}
	
	val=simps_sum(sample_lap,2*n+1,step);
	//printf("%.3e\n",val);
	return(val);
}

double grad_k_2(double k,double step){
	double val=0;
	double r=0;
	double kr=0;
	for(int j=0;j<(2*n+1);j++){
		r=R_MIN+j*step;
		kr=k*(r);
		val=( dbesj0_(&kr)/r - dbesj1_(&kr) *k) * sample[j] ;
		sample_lap[j]=val;	
	}
	
	val=simps_sum(sample_lap,2*n+1,step);
	return(val);
}

//////////////////////////////////////////////////////////////////////////
double saturation(double step,double* sudpar,double Q2){
	double k_step=0.1;
	double k_min=0.4;
	double k=k_min;
	double prev=1;
	double val;
	
	for(int i=0;i<6;i++){
		//printf("%.3e\t%.3e\n",k_min,k_step );
		k=k_min;
		
		for(int j=0;j<100;j++){
			val=grad_k(k,step,sudpar,Q2);
			//printf("%d : %.3e, %.3e, %.3e, %.3e\n",j ,k,k_min, val, prev );
			if(j!=0 ){
				if(prev*val<0){
					k_min=k-k_step;
					//printf("%d : %.3e, %.3e, %.3e, %.3e\n",i ,k,k_min, val, prev );
					break;
				}
			}
			if(fabs(val)>fabs(prev)){
				//printf("error\t%.3e\t%.3e\n",val,prev);
			}
			//printf("%f, %f, %f\n",k, val, prev );
			prev=val;
			k+=k_step;
		}
		if(fabs(prev)<1.0e-5){
			break;
		}
		if(i!=5){
			k_step*=0.05;
		}
		
	}
	//printf("\nk=%.3e, val= %.3e, prev= %.3e\n\n",k, val, prev );
	
	k=k - k_step*fabs(val/(prev-val) );//weighted mid point 	
	return(k);
}
