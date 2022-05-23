//#include<stdio.h>
//#include<math.h>

//#include"./dipole-cross-section.h"
//#include"./photon-wave-function.h"
//#include"./simpson-integral.h"
//#include"./constants.h"


//#define FLAVOUR 2 

double sigma_integrand(double R,double z,double * par[2]){

	//R is r/(1+r) such that integration between (0, inf) -> (0,1).
	double r=R/(1-R);
	double jacob=pow(1-R,-2);	
	double x=(*(*par));
	double Q2=(*( (*par)+1 ));

       	double* param=( *(par+1) );
	double value=0;
	

	#if FLAVOUR==0
       		value=sigma(r,x,Q2,param,MODEL,'l')*psisq_f(r,z,Q2,  'l' );
                value+=sigma(r,x,Q2,param,MODEL,'s')*psisq_f(r,z,Q2, 's' );
	#else 
		#if FLAVOUR==1
		       	value=sigma(r,x,Q2,param,MODEL,'l')*psisq_f(r,z,Q2,  'l' );
                	value+=sigma(r,x,Q2,param,MODEL,'s')*psisq_f(r,z,Q2, 's' );
                	value+=sigma(r,x,Q2,param,MODEL,'c')*psisq_f(r,z,Q2, 'c' );
		#else 
			#if FLAVOUR==2
        		value=sigma(r,x,Q2,param,MODEL,'l')*psisq_f(r,z,Q2,  'l' );
			value+=sigma(r,x,Q2,param,MODEL,'s')*psisq_f(r,z,Q2, 's' );
			value+=sigma(r,x,Q2,param,MODEL,'c')*psisq_f(r,z,Q2, 'c' );
			value+=sigma(r,x,Q2,param,MODEL,'b')*psisq_f(r,z,Q2, 'b' );
			#endif
		#endif
	#endif
	return(jacob*r*value);
}





double sigma_DIS(double x,double q2,double y, double * par) {

	double sigma_sum;

	double var[2]={x,q2};
	double *param[2]={var,par};
	double res=0;

simpson2d(& sigma_integrand, param, 1.0e-10,1.0-1.0e-10,0.0,1.0,&res);
	
	return res;
}	

