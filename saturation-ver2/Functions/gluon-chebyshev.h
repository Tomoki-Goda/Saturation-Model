//#include<math.h>
//#include<stdlib.h>
//#include<stdio.h>


//#include"./control_tmp.h"
//#include"control.h"
//#include"./control-default.h"
//#include"./constants.h"

#include"./gluons.h"
#include"./chebyshev-1.h"

#define X_DEGREE 25 
#define Q2_DEGREE 27
 
#ifndef SKIP_SAME
	#define SKIP_SAME 0
#endif

static int ERROR_FLAG = 1;
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
static double COEFF[X_DEGREE*Q2_DEGREE];
static const unsigned DIM=2;
static const unsigned DEGREE[2]={X_DEGREE,Q2_DEGREE};

//static const double X_LIM[2]={1.0e-8,1} ;//x upper lim is not 0.01 for it is x_mod... that can even be >1 for small Q/m but in practice 2m/Q is bound by pair production threshold.
static const double X_LIM[2]={1.0e-8,1.25} ;//x upper lim is not 0.01 for it is x_mod... that can even be >1 for small Q/m but in practice 2m/Q is bound by pair production threshold.
static const double Q2_LIM[2]={LQCD2*3 , 5.0/(R_MIN*R_MIN) };

/////////////////////////////////////////

///////////////  format  //////////////////
double eval_xg(const double *var, const double *par ){
	set_xg_parameter( par[0] , par[1]);
	//double x= change_var_revert(X_LIM[0],X_LIM[1],var[0]);
	double x= change_var_revert_log(X_LIM[0],X_LIM[1],var[0]);
	double Q2= change_var_revert_log(Q2_LIM[0],Q2_LIM[1] ,var[1]);
	
	double val=alpha_s(Q2 )* xgpdf(x,Q2);
	return( val);
}

////////////// approximation /////////////
int comp( double a, double b ,double tolarence){
	return( (int)(fabs(a-b)/(( fabs(a)+fabs(b) )*tolarence)) );
}
void approx_xg(const double * par){
	//this par is sigmapar+1. shoud be written so in read-and-fit.
	
//	if(par[3]/(R_MIN*R_MIN) >  Q2_LIM[1]  ){
//		printf("Q2_LIM too small: C=%.3e r_min=%.3e Q2_LIM=%.3e\t " ,par[3],R_MIN, Q2_LIM[1] );
//		Q2_LIM[1]=1.1*par[3]/(R_MIN*R_MIN);
//		printf("New Q2_LIM=%.3e\n",Q2_LIM[1] );
//		goto re_evaluate;
//	}
	//printf("%.3e %.3e %.3e",par[2],par[3],par[4] );
	
//	Q2_LIM[0]=par[2]/(1.1*par[3]*par[3]);
//	Q2_LIM[1]=1.1*par[2]/(R_MIN*R_MIN);
	//printf("Q2 for xg() between %.3e %.3e\n", Q2_LIM[0] , Q2_LIM[1]  );
//	if(Q2_LIM[0]<LQCD2){
//		printf("Q2 might get less than Lambda QCD: %.3e\n",Q2_LIM[0]);
//	}
	
	//static double param[2];
	
	//clock_t time_xg=clock();
	//printf("*************************       xg(x,Q2)          **************************\n");
//#if SKIP_SAME==1
//	if((comp(*(par+1),param[1] , 1.0e-10)==0)){
//		if((comp(*(par),param[0] , 1.0e-10)==0)){
//			//printf("skip\n");				//both A_g and lambda_g are not changed.
//		}else{	
//			for(unsigned i=0;i<(X_DEGREE*Q2_DEGREE); i++){
//				*(coeff+i)=((*(par))/param[0])*(*(coeff+i));
//			}
//			//printf("skip - normalize\n");			//only A_g is different
//			*(param)   = *(par);
//		}
//		
//	}else{
//		*(param)   = *(par);
//		*(param+1) = *(par+1);
//#endif
//		re_evaluate:
		cheb_coeff( &eval_xg, par,DEGREE,DIM, COEFF  );	// both A_g and lambda_g are different. 
//#if SKIP_SAME==1
//	}
//#endif
	//time_xg-=clock();
	//printf(  "*********************   Chebyshev  done, %.3e seconds    **************************\n", -((double)time_xg)/CLOCKS_PER_SEC);
	ERROR_FLAG=0;//indicate approxmation is done,
	
}

///////////// result /////////////////
double xg_chebyshev(double  x,double q2){
	if(ERROR_FLAG!=0 ){
		printf("approx_xg has to be run first!!!");
		getchar();
		//char ch;
		//scanf("%c",&ch);
	}
	if(q2<Q2_LIM[0]){
		//printf(" q2 too small %.2e\n",q2);
		return 0;
		//q2=Q2_LIM[0];
	}else if(q2>Q2_LIM[1]){
		//printf(" q2 too large %.2e\n",q2);
		return 0;
		//q2=Q2_LIM[1];
	}
	
	if(x<X_LIM[0]){
		//printf(" x too small %.2e\n",x);
		return 0;
		//x=X_LIM[0];
	}else if(x>X_LIM[1]){
		//printf(" x too large %.2e\n",x);
		return 0;
		//x=X_LIM[1];
	}
	
	
	double args[2];
	//(*args)=change_var_compactify(X_LIM[0],X_LIM[1],x);
	(*args)=change_var_compactify_log(X_LIM[0],X_LIM[1],x);
	(*(args+1))=change_var_compactify_log(Q2_LIM[0],Q2_LIM[1],q2);
	
	double val = chebyshev(DEGREE, DIM, COEFF, args );
	
	
	//char ch;
	//printf("value = %f x = %f  q2 = %f\n",val,x,q2);
	if(val<0){
		//printf(" alpha x g(%.2e,%.2e)=%.3e\n",x,q2,val);
		val=0;
	}
	//printf("%f arg1 = %f  arg2 = %f\n",val, *args,*(args+1));
	//scanf("%c",&ch);
	return val;
} 
