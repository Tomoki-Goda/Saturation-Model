//#include<math.h>
//#include<stdlib.h>
//#include<stdio.h>
//#include"./constants.h"
//#include"./control-default.h"
#include"./gluons.h"
#include"./chebyshev-1.h"

#define X_DEGREE 50
#define Q2_DEGREE 50
 
#ifndef SKIP_SAME
	#define SKIP_SAME 1
#endif

///////just to clarify/////////
//extern void cheb_coeff(double func(double * vec,double* par) , double * par,const unsigned *degree,unsigned dim, double* coeff  );
//extern double chebyshev(const unsigned *degree,unsigned dim, double* coeff ,double* args );

//extern double change_var_revert(double min,double max, double val)
//extern double change_var_compactify(double min,double max, double val)
////////////
//extern void set_xg_parameter(double A_g,double lambda_g);
/////////////////////////////
extern double alpha_s(double mu2 );

////////////////Global///////////////////
static double coeff[X_DEGREE*Q2_DEGREE];
static  const unsigned dim=2;
static  unsigned degree[2]={X_DEGREE,Q2_DEGREE};

static double xlim[2]={1.0e-7,1.0} ;//x upper lim is not 0.01 for it is x_mod... 
static double q2lim[2]={0.5,1.0e+10};
/////////////////////////////////////////

///////////////  format  //////////////////
double eval_xg( double *var,double *par ){
	set_xg_parameter( par[0] , par[1]);
	double x= change_var_revert_log(xlim[0],xlim[1],var[0]);
	double Q2= change_var_revert_log(q2lim[0],q2lim[1] ,var[1]);
	//double x=(var[0]+1)/2;
	//double Q2=(1+var[1])/(1-var[1]);
	double val=alpha_s(Q2 )* xgpdf(x,Q2);
	return( val);
}

////////////// approximation /////////////
int comp( double a, double b ,double tolarence){
	return( (int)(fabs(a-b)/(( fabs(a)+fabs(b) )*tolarence)) );
}
void approx_xg(double * par){
	static double param[2];
	clock_t time_xg=clock();
	//printf("*************************       xg(x,Q2)          **************************\n");
#if SKIP_SAME==1
	if((comp(*(par+1),param[1] , 1.0e-10)==0)){
		if((comp(*(par),param[0] , 1.0e-10)==0)){
			printf("skip\n");				//both A_g and lambda_g are not changed.
		}else{	
			for(unsigned i=0;i<(X_DEGREE*Q2_DEGREE); i++){
				*(coeff+i)=((*(par))/param[0])*(*(coeff+i));
			}
			printf("skip - normalize\n");			//only A_g is different
			*(param)   = *(par);
		}
		
	}else{
#endif
		*(param)   = *(par);
		*(param+1) = *(par+1);
		cheb_coeff( &eval_xg, par,degree,dim, coeff  );	// both A_g and lambda_g are different. 
#if SKIP_SAME==1
	}
#endif
	time_xg-=clock();
	printf(  "*********************   Chebyshev  done, %.3e seconds    **************************\n", -((double)time_xg)/CLOCKS_PER_SEC);
	
}

///////////// result /////////////////
double xg_chebyshev(double  x,double q2){
	double args[2];
	(*args)=change_var_compactify_log(xlim[0],xlim[1],x);
	(*(args+1))=change_var_compactify_log(q2lim[0],q2lim[1],q2);
	
	double val = chebyshev(degree, dim, coeff, args );
	
	
	char ch;
	//printf("%f x = %f  q2 = %f\n",val,x,q2);
	if(val<0){
		//printf(" alpha x g(%.2e,%.2e)=%.3e\n",x,q2,val);
		val=0;
	}
	//printf("%f arg1 = %f  arg2 = %f\n",val, *args,*(args+1));
	//scanf("%c",&ch);
	return val;
} 
