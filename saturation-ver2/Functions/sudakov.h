//#include<stdio.h>
//#include<math.h>


//#include"./control_tmp.h"
//#include"control.h"
//#include"./control-default.h"
//#include"constants.h"

//#include"dipole-cross-section.h"

#include"./mu2.h"


//extern double sudakov(double ,double, double* );// in dipole-cross-section
extern double sigma_gbw(double, double, double, double*);
extern double sigma_bgk(double, double, double, double*);

///////////////////////////////////////////////////////////////////////////////////////////////
extern double dgauss_(double(* )(double*),double*,double*,double*); 
extern double dgquad_(double(* )(double*),double*,double*,int*);
///////////////////////////////////////////////////////////////////////////////////


double sudakov(double r, double mu2,double* par) {
	//pertubative+non-perturbative sudakov
	double C= (*(par));
	double r_max=( *(par+1));
	double g1=(*(par+2));
	double g2=(*(par+3));
	
	double mub2=rmu2(r ,par );

	double val=0.0;

#if SUDAKOV>=2	
	val+=(g1/2) * pow(r,2) + (g2/4) * log(1+pow(r/r_max,2) )*log( mu2 /pow(Q0 ,2) );	
#endif
	
#if THETA_OFF==0
	if (mu2 < mub2) {//ensures that lower limit of integral is smaller than upper limit...
		return(val);
	}
#endif

#if SUDAKOV>=1
	if (mub2 < LQCD2){
		printf("\nsudakov:: mu_b is too low!!!\n%f\t%f\n\n",C,r_max);
	}
	if (mu2 < LQCD2 ) {
		//treatment here is unsure, this is dependent on how we treat the case Q2<mu2b
		//return(val);
		//printf("\nsudakow:: warning non perturbative region: Q2=%f\n\n",mu2);
		return(0.0);
	}
	
	double b0 = (11*CA-2*NF)/12;
	double L_mu_l=log(mu2/LQCD2);
	double L_mub_l=log(mub2/LQCD2);
	
	val += CA/(2*b0*PI)*(L_mu_l*log(L_mu_l /L_mub_l)-(L_mu_l-L_mub_l ));
	//double val = CA/(2*b0*PI)*(log(mu2/LQCD2)*log(log(mu2/LQCD2)/log(mub2/LQCD2))-log(mu2/mub2));
#endif	

    return val;
}


//////////////////////////////////jacobian ////////////////////////////////////////////
double dsdr(double r,double Q2, double* sudpar ){
#if SUDAKOV==0
	return 0.0;
#endif
	
	double mu2=rmu2(r,sudpar);
	double jacfirst=rmu2_jac_first(r,sudpar);
	
	if(mu2>Q2){ 
		return(0.0);
	}
	if(mu2<LQCD2){
		printf( "drds:: mu2 out of range");
	}
	double b0= (11.0*CA-2.0*NF)/12.0;
	double alpha =1.0/(b0* log( mu2/LQCD2) );
	double val= - jacfirst* alpha*log(Q2/mu2)/mu2;
	return(((double)CA)/(2.0*PI) * val);
}

double ddsdrdr(double r,double Q2,double * sudpar){
#if SUDAKOV==0
	return 0.0;
#endif
	double mu2=rmu2(r,sudpar);
	double jacfirst=rmu2_jac_first(r,sudpar);
	double jacsecond=rmu2_jac_second(r,sudpar);
	
	double b0 = (11.0*CA-2.0*NF)/12.0;
	if(mu2>Q2){
		return(0.0);
	}
	if(mu2<LQCD2){
                printf( "drds:: mu2 out of range");
        }

	double val;

	double logQmu=log(Q2/mu2);
	double alpha =1.0/(b0* log( mu2/LQCD2) );

	//derivative s wrt mu2 with jacobians
	double val1	= - jacsecond * alpha*logQmu/mu2;	
	double val2	=   jacfirst*jacfirst * alpha/pow(mu2,2) * ( 1 + logQmu*(1+b0*alpha));
	
	val=val1+val2;
	return(CA/(2*PI) * val); 
}

////////////////////////////Integration/////////////////////////////
static double VAR[3];
static double *SIGPAR;
static double *SUDPAR;


double integrand( double * r_ptr){
	double val;
	double r=*r_ptr;
	double  R=VAR[0];
	double  x=VAR[1];
	double  Q2=VAR[2];
	
	//double C=PAR[3];
	//double rmax=PAR[4];
	//double C2sq=sqrt(PAR[5]);
	
	double mu2=rmu2(r, SUDPAR);
	
	if(Q2<mu2){
		return 0.0;
	}
	if(mu2<LQCD2){
		printf("integrand::mu2 out of range");
	}
		
	double dels=dsdr(r, Q2, SUDPAR);
	double deldels=ddsdrdr(r,Q2, SUDPAR);

	double sigma_0=*(SIGPAR);
	double sud=sudakov(r,Q2,  SUDPAR );
	
	double logrr=log(R/r);
	
	//printf("%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t" ,dels,deldels,sigma_0,sud,logrr);
	
	val=(2.0-logrr)*dels+r*logrr*(dels*dels-deldels) ;
	//printf("%.2e\t" ,val);
	val*=exp(-sud)*BASE_SIGMA(r,x,Q2,SIGPAR);
	//printf("%.2e\n" ,val);
	
	return val;
}




double integral_term(double r, double x, double Q2, double * sigmapar, double* sudpar){
	//clock_t time=clock();
	
	VAR[0]=r;
	VAR[1]=x;
	VAR[2]=Q2;
	SIGPAR=sigmapar;
	SUDPAR=sudpar;

	double result=0;
	double rmin=1.0e-5;
	//double C=*(par+3);
	//double rmax=*(par+4);
	//double C2sq=sqrt(PAR[5]);
	double mu2=rmu2(r, SUDPAR);
	
	//printf("1/rmin^2 =%.1e 1/r^2=%.1e \n" ,invrmin2,1.0/(r*r));
	if( mu2>Q2 ){ 
		return(0.0);
	}
	rmin= pow( rmin2(Q2, SUDPAR ) ,0.5);
	//printf("rmin =%.1e r=%.1e \n" ,rmin,r);
	double N=DGAUSS_PREC;
	//int N=96;
	//rmin=1.0/pow(invrmin2,0.5);
	//printf("rmin= %.1e\n",rmin);
	//result=dgquad_(&integrand,&rmin,VAR,&N);
	result=dgauss_(&integrand,&rmin,VAR,&N);	
	//printf("%f\n",result);
	//time-=clock();
	//printf("%.2e seconds\n", -((double)time)/CLOCKS_PER_SEC);
	//for (unsigned i=0;i<5;i++){
	//	printf("%.2e\t",sigmapar[i]);
	//}
	//
	//for (unsigned i=0;i<2;i++){
	//	printf("%.2e\t",sudpar[i]);
	//}
	//printf("result %f\n",result);
	
	return result;
	
}

double sigma_s(double r, double x, double Q2, double * sigmapar, double* sudpar){
	//for (unsigned i=0;i<5;i++){
	//	printf("%.2e\t",sigmapar[i]);
	//}
	//printf("\n");
	//for (unsigned i=0;i<2;i++){
	//	printf("%.2e\t",sudpar[i]);
	//}
	//printf("\n");
	
	double sud=exp(-sudakov(r,Q2, sudpar ) );
	
	double val=sud*BASE_SIGMA(r,x,Q2, sigmapar );
	//printf("val sud val %.2e \t %.2e \t",val,sud);
	val+=integral_term(r,x,Q2,sigmapar,sudpar);
	//printf("%.2e \n",val);
	return val;
}





