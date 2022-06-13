#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include"./constants.h"
#include"./control-default.h"
#include"./gluons.h"
#include"./chebyshev-1.h"

#define X_DEGREE 20
#define Q2_DEGREE 20

///////just to clarify/////////
//extern void cheb_coeff(double func(double * vec,double* par) , double * par,const unsigned *degree,unsigned dim, double* coeff  );
//extern double chebyshev(const unsigned *degree,unsigned dim, double* coeff ,double* args );
////////////
//extern void set_xg_parameter(double A_g,double lambda_g);
/////////////////////////////


////////////////Global///////////////////
static double coeff[X_DEGREE*Q2_DEGREE];
static const unsigned dim=2;
static const unsigned degree[2]={X_DEGREE,Q2_DEGREE};
/////////////////////////////////////////

//double change_var(double min,double max, double val){
//	if((val>1)||(val<-1)){
//		printf("wring input for change_var");
//		return 0.0;
//	}
//	return max*pow((min/max),(1.0-val)/2);
	//max*(min/max)^((1-x)/2) =[min,max] for x=[-1,1]	
//}
///////////////////////////Change of variables to -1 1 ///////////////////////
//double change_var_compactify(double min,double max, double val){
double change_var_revert(double min,double max, double val){
	//for val =[-1,1]  return value between min and max
	if((val>1)||(val<-1)|| (min>max)){
		char str[100];
		printf("wrong input for change_var_revert\n val=%f\t [%f, %f] \n",val, min,max);
		printf("%f\t%f\t %f \n type something and enter to continue\n",val, min,max);
		scanf("%s ",str);
		printf("continue");
		//return 0.0;
	}
	return ( (val/2) *(max-min)+ (max+min)/2 );
}
//double change_var_revert(double min,double max, double val){
double change_var_compactify(double min,double max, double val){
	//for val =[min,max]  return value between -1 and 1
	//if((val>max)||(val<min)|| (min>max)){
		//char str[100];
		//printf("wrong input for change_var_compactify\n val=%f\t [%f, %f] \n",val, min,max);
		//printf("%f\t%f\t %f \n type something and enter to continue\n",val, min,max);
		//scanf("%s ",str);
		//printf("continue");
		//return 0.0;
	//}
	return (2*((val-min)/(max-min)) -1);
	
}

///////////////  format  //////////////////
double eval_xg( double *var,double *par ){
	set_xg_parameter( par[0] , par[1]);
	double x= change_var_revert(1.0e-7,1.0,var[0]);
	double Q2= change_var_revert(1.0e-7,1000 ,var[1]);
	//double Q2= change_var_revert(1.0e-7,1-(1.0e-7) ,var[1]);
	//Q2=Q2/(1-Q2);
	
	//double Q2;= change_var_revert(1.0e-2,1.0e+3,var[1]);
	//printf("eval xg\n");
	//char str[10];
	//scanf("%s",str );
	double val=xgpdf(x,Q2);// var has to be put in the right region ???
	return( val);
}



////////////// approximation /////////////
void approx_xg(double * par){
	printf("\n*************************************** xg(x,Q2) **************************\n\n");
	cheb_coeff( &eval_xg, par,degree,dim, coeff  );
	printf(  "***************************************   Done   **************************\n");
}

///////////// result /////////////////
double xg_chebyshev(double  x,double q2){
	double args[2];
	
	
	(*args)=change_var_compactify(1.0e-7,1.0,x);
	
	(*(args+1))=change_var_compactify(1.0e-7,1000,q2);
	//(*(args+1))=change_var_compactify(1.0e-7,1-(1.0e-7),q2/(1+q2) );
	
	char ch;
	double val = chebyshev(degree, dim, coeff, args );
	//printf("%f x = %f  q2 = %f\n",val,x,q2);
	//printf("%f arg1 = %f  arg2 = %f\n",val, *args,*(args+1));
	//scanf("%c",&ch);
	return val;
} 
