
//extern double Kahn_list_sum(const double *arr,int len);
//#include"./kahnsum.h"
//#include"./critical-line.h"
//#include"../Functions/sudakov.h"
#include"../Functions/Kahn.h"
//#include"./main.h"
#include"clenshaw.h"
extern double dbesj0_(double*);
extern double dbesj1_(double*);
extern double mod_x(double, double,int);
extern double SIGMA(double , double ,double , const double *, const double*);
extern double BASE_SIGMA(double , double ,double , const double *);
extern double laplacian_sigma(double ,double,double,double*,double *);

extern int parameter(const double*,double*,double*);
extern void approx_xg(const double *);


static const int n=750;
static double sample[2*750+1 +2]={0};
static double sample_lap[2*750+1]={0};
static double summand[2*750+1]={0};

static struct parameters{
double x;
double Q2;
double *sigpar;
double k;
double *sudpar;
double *mu2;
} PARAM;

//static double lap_for_int(double*r){
//	return(laplacian_sigma(PARAM.x, *r, PARAM.Q2,PARAM.sigpar));
//}
static double integrand_af(double* r,void* ){
	double r2=(*r)/(1-*r);
	double jac=1.0/pow(1-*r,2);

	double kr=PARAM.k * r2;
	double val=jac * 2*PI*dbesj0_(&kr)*(r2)*laplacian_sigma(r2,PARAM.x, PARAM.Q2,PARAM.sigpar,PARAM.sudpar);

 	//printf("%.5e\tr=%.5e\n",val,*r);
	return val;
}

extern double dgauss_(double (*)(double*),double*,double*,double*);
extern double dgquad_(double(*)(double *), double*, double*,int*);
extern double dadapt_(double(*)(double *), double*, double*,int*,double*,double*,double*,double*);
//extern double dclenshaw(double (*)(double *), double ,double, double);

double af(double x,double k,double q2,double * sigpar,double *sudpar){
	PARAM.x=x;
	PARAM.Q2=q2;
	PARAM.sigpar=sigpar;
	PARAM.sudpar=sudpar;
	PARAM.k=k;
	double eps=1.0e-8;
	double epsabs=0;
	int seg=1;
	double val, err;
	double lim[2];
	//lim[0]=R_MIN;
	//lim[1]=R_MAX;
	lim[0]=R_MIN/(1+R_MIN);
	lim[1]=((double)R_MAX)/(1+R_MAX);
	double* dummy;	
#if IIM==1
	val=0;
	double Qs=pow(1.632e-5/x,0.2197/2);
	double r_bound=2.0/Qs;
	r_bound=r_bound/(1+r_bound);
	printf("%.3e %.3e %.3e\n", lim[0],r_bound, lim[1]);
	if(r_bound<lim[0]){
		r_bound=lim[0];
		printf("skip\n " );
	}else{
		val+=dclenshaw(&integrand_af,(void*)dummy,lim[0],r_bound-1.0e-3,eps);
		printf("val1 %.3e\n" ,val);
	}
	if(r_bound>lim[1]){
		r_bound=lim[1];
		printf("skip\n " );
	}else{
		val+=dclenshaw(&integrand_af,(void*)dummy,r_bound+1.0e-3,lim[1],eps);
		printf("val2 %.3e\n" ,val);
	}
#else
	val=dclenshaw(&integrand_af,(void*)dummy,lim[0],lim[1],eps);
#endif
	//dadapt_(&integrand_af,lim,lim+1,&seg,&eps,&epsabs,&val,&err);
	//double val=dgquad_(&integrand_af,lim,lim+1,&N);
	printf("af= %.5e+-%.5e \tk^2=%.5e \n",val,err,k*k);
	return val;
}


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
//#if ((MODEL==0||MODEL==1)||(PHI==1))
//		val=SIGMA(r,x,Q2,sigpar,sudpar);
//#else 
		val=BASE_SIGMA(r,x,Q2,sigpar);
//#endif
			
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
	val=Kahn_list_sum(summand,2*n+1);
	val*=(step/3);
	return val;
}

/////////////////////////////////////////////////////////////////////////////////
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
			
			//val*=exp_sud(r,mu2_arr[0],q2);
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
			//signal=compute_mu2(r, sudpar, mu2_arr,1);//compute mu2

			val=(sample[j]-2*sample[j-1]+sample[j-2])/(step*step) + (sample[j]-sample[j-2])/(2*step* (r-step));
			//val*=exp_sud(r,mu2_arr[0],q2);
			//val=sample[2*n+1]- sample[j-1];
			kr=k*(r-step);
			
			val*=((r-step)*(r-step)*dbesj1_(&kr));//This line is the only difference
			
			sample_lap[j-2]=val;	

		}
	}
	
	val=simps_sum(sample_lap,2*n+1,step);
	return(val);
}


//////////////////////////////////////////////////////////////////////////
double saturation(double step,double* sudpar,double Q2){
	double k_step=0.1;
	double k_min=0.3;
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

