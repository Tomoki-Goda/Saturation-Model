#include<math.h>
#include<stdlib.h>
#include<stdio.h>

#include"./kahnsum.h"
//extern double k_group_sum(double *, int);
double boole_coeff(){
	static int index=0;
	static double coeff_arr[8]={32,12,14,32,32,14,12,32};
	double val=coeff_arr[index];
	(index==7)?(index=0):(index++);
	return( val);
}

//double ARR[10000];

double boole_update(double (*func)(double, double**), double** par, double prev,double min, double *step, int *len){
	/////////////////////////////////////
	//len is N-1 for N point integration.　	multiple of 4.
	//step=(max-min)/N
	////////////////////////////////////
	double val=0;
	(*step)/=3;
	double arr[2*(*len)];

	for(int i=0;i<(*len);i++){
		arr[2*i]=boole_coeff()*(*func)(min+((3*i+1) *(*step)),par);
		arr[2*i+1]=boole_coeff()*(*func)(min+((3*i+2) *(*step)),par);
		//val+=boole_coeff()*(*func)(min+((3*i+1) *(*step)),par);
		//val+=boole_coeff()*(*func)(min+((3*i+2) *(*step)),par);
	}	
	val=k_group_sum(arr,2*(*len));
	(*len)*=3;
	return(prev+val);
}

double boole_integral(
		double (*func)(double, double**),double **par, 
		double min, double max,
		double eps_abs,double eps_rel, 
		int min_N, int max_N,double *res,double *error){
	if(4*pow(3,max_N)>1.0e+5){
		printf("%.0f is too large Reduce max_N \n",pow(3,max_N)*4+1 );
	}
	double step=(max-min)/4;
	double factor=(2*step)/45;
	//////////////initial interation/////////////////
	double val=7*((*func)(min,par)+(*func)(max,par)) 
		+ 32*( (*func)(min+step,par)+(*func)(min+3*step,par)) 
		+ 12*(*func)(min+2*step,par);
	////////////////////////////////////////////////
	//printf("%.5e\n",factor*val);	
	int len=4;
	double new,error_abs;
	for(int i=0;i<max_N;i++){
		factor/=3;
		new=boole_update(func,par,val,min,&step,&len);
		error_abs=factor*fabs(val-(new/3));
		val=new;
		if((i>=(min_N+1)) &&((error_abs<eps_abs)||(fabs(error_abs/(factor*val))<eps_rel)  ) ){
			break;
		}
		//printf("%i\t %.5e +- %.5e  %.5e %.5e \n",i, factor*val,error_abs,fabs(error_abs/(factor*val)) , 4*pow(3,(i+1))+1  );		
	}
	*res=step*(2*val)/45;
	*error=error_abs;
	//if(fabs(error_abs/(factor*val))>eps_rel){
	//	printf("%.5e +- %.5e  %.5e \n", factor*val,error_abs,fabs(error_abs/(factor*val)) );		
	//}
	return(*res);
}



///////////Test////////
/*double integrand(double r, double** par){
	return(pow(r,2)*exp(-r*r*r*r));
}
int main(){
	double **par;
	double val,err;
	boole_integral(&integrand,par,0,50,0,1.0e-15,3,7,&val,&err);
	printf("%.5e +- %.5e\n",val,err);
	return(0);
}

*/











