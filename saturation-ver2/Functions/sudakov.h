//#include<stdio.h>
//#include<math.h>


//#include"./control_tmp.h"
//#include"control.h"
//#include"./control-default.h"
//#include"constants.h"

//#include"dipole-cross-section.h"

#include"./mu2.h"


//extern double sudakov(double ,double, double* );// in dipole-cross-section
//extern double BASE_SIGMA(double, double, double, double*);
extern double sigma_gbw(double, double, double, const double*);
extern double sigma_bgk(double, double, double, const double*);

///////////////////////////////////////////////////////////////////////////////////////////////
extern double dgauss_(double(* )(double*),double*,double*,double*); 
extern double dgquad_(double(* )(double*),double*,double*,int*);
extern double dadapt_(double(* )(double*),double*,double*,int*,double*,double* ,double*,double*);

static double VAR[3];
static const double *SIGPAR;
static const double *SUDPAR;

///////////////////////////////////////////////////////////////////////////////////

double sudakov_p(double mub2, double Q2,const double* par) {
		
	if (Q2 < mub2) {//ensures that lower limit of integral is smaller than upper limit...
		return(0.0);
	}
	if (mub2 < LQCD2){
		printf("\nsudakov:: mu_b is too low!!!\t mub2= %.3e\n",mub2);//\n%f\t%f\n\n",C,r_max);
		getchar();
	}
	if (Q2 < LQCD2 ) {
		printf("\nsudakov:: Q2 is too low!!!\t Q2= %.3e\t mub2= %.3e\n",Q2,mub2);
		//printf("%.3e\t%.3e\t%.3e\n",par[0],pow(r,-2.0),pow(par[1],-2.0));
		getchar();
	}
	double b0 = (11*CA-2*NF)/12;
	double L_Q_l=log(Q2/LQCD2);
	double L_mub_l=log(mub2/LQCD2);
	
	//double val=pow( Q2/mub2 * pow (L_Q_l/L_mub_l,-L_Q_l  ), CA/(2.0*PI *b0) );
	double val=CA/(2*b0*PI)*(L_Q_l*log(L_Q_l /L_mub_l)-(L_Q_l-L_mub_l ));
	//val=exp(-val);
	return val;
}
double sudakov_np(double r,double mub2, double Q2,const double* par) {
	double C= (*(par));
	//double r_max=( *(par+1));
	double g1=(*(par+2));
	double g2=(*(par+3));
	
	//double mub2=rmu2(r ,par );
	if (mub2 < LQCD2){
		printf("\nsudakov:: mu_b is too low!!!\t mub2= %.3e\n",mub2);
	}
	
	double val=g1/2 * r*r;
	
	//if (Q2 < Q0 ) {
	//	return(val);
	//}
	
	val+=(g2/4) * log(mub2*(r*r)/C )*log( Q2 /pow(Q0 ,2) );
	
	return val;
}
	
//////////////////////////////////jacobian ////////////////////////////////////////////

double dsnpdr(double r, const double * mub2_arr, double Q2,const double* sudpar){
#if SUDAKOV<=1 
	return 0.0;
#endif 	
	double rmax=sudpar[1];
	double g1=sudpar[2];
	double g2=sudpar[3];
	double val;
	
	double mu2=mub2_arr[0];
	double jacfirst=mub2_arr[1];//rmu2_jac_first(r,sudpar);
	
	val=g1*r + (g2/4)* (2.0/r+jacfirst/mu2)*log(Q2/(Q0*Q0)); 
	
	return(val);

}


double ddsnpdrdr(double r,const double * mub2_arr, double Q2,const double* sudpar){
#if SUDAKOV<=1 
	return 0.0;
#endif 	
	double rmax=sudpar[1];
	double g1=sudpar[2];
	double g2=sudpar[3];
	double val;
	double mu2=mub2_arr[0];
	double jacfirst=mub2_arr[1];
	double jacsecond=mub2_arr[2];
	
	val=g1+(g2/4)*(- 2.0/(r*r) +jacsecond/mu2 -(jacfirst*jacfirst)/(mu2*mu2))*log(Q2/(Q0*Q0)) ; 

	return(val);

}



double dsdr(const double * mub2_arr,double Q2, const  double* sudpar ){
#if SUDAKOV==0
	return 0.0;
#endif
	
	double mu2=mub2_arr[0];
	double jacfirst=mub2_arr[1];

	if(mu2>Q2){ 
		return(0.0);
	}
	if(mu2<LQCD2){
		printf( "drds:: mu2 out of range");
	}
	double b0= (11.0*CA-2.0*NF)/12.0;
	double alpha =1.0/(b0* log( mu2/LQCD2) );
	
	double val= - jacfirst* alpha*log(Q2/mu2)/mu2;

	val*=((double)CA)/(2.0*PI);

	return(val);
}

double ddsdrdr(const double * mub2_arr,double Q2,const double * sudpar){
#if SUDAKOV==0
	return 0.0;
#endif
	double mu2=mub2_arr[0];
	double jacfirst=mub2_arr[1];
	double jacsecond=mub2_arr[2];
	//double mu2=rmu2(r,sudpar);
	//double jacfirst=rmu2_jac_first(r,sudpar);
	//double jacsecond=rmu2_jac_second(r,sudpar);
	
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
	val*= ((double)CA)/(2*PI) ; 

//#if SUDAKOV==2
//	val+=ddsnpdrdr(r,Q2,sudpar);
//#endif
	return(val);

}
////////////////////////////Integration/////////////////////////////


double integrand( double * r_ptr){
	double val;
	double r=*r_ptr;
	double  R=VAR[0];
	double  x=VAR[1];
	double  Q2=VAR[2];
	
	//double C=PAR[3];
	//double rmax=PAR[4];
	//double C2sq=sqrt(PAR[5]);
	double dels,deldels;
	
	//double mu2=rmu2(r, SUDPAR);
	double mu2_arr[3];
	compute_mu2(r, SUDPAR, mu2_arr,3);//compute mu2, dmu2/dr dd mu2/drdr
	
	if((mu2_arr[0]) <LQCD2){
		printf("integrand::mu2 out of range\n");
		printf("R %f  Q2 %f x %f r %f mu2 %f\n" , R, Q2, x, r , (mu2_arr[0]) );
		getchar();
	}
	
	if((mu2_arr[0]) >Q2){
		dels=dsdr(mu2_arr, Q2, SUDPAR);
		deldels=ddsdrdr(mu2_arr,Q2, SUDPAR);
#if SUDAKOV==2
		dels+=dsnpdr(r,mu2_arr, Q2, SUDPAR);
		deldels+=ddsnpdrdr(r,mu2_arr,Q2, SUDPAR);
#endif
	}else{
#if SUDAKOV==2
		dels=dsnpdr(r,mu2_arr, Q2, SUDPAR);
		deldels=ddsnpdrdr(r,mu2_arr,Q2, SUDPAR);
#else
		return(0.0);
#endif
	}
	//if(Q2<mu2){
	//	printf("check lower limit of r ");
	//	getchar();
	//	return 0.0;
	//}
	double sud=0.0;
	
#if SUDAKOV>=1
	sud+=sudakov_p(mu2_arr[0] ,Q2, SUDPAR);
#endif
#if SUDAKOV==2
	sud+=sudakov_np(r, mu2_arr[0],Q2, SUDPAR);
#endif
	double logrr=log(R/r);
	
	//printf("%.2e\t%.2e\t%.2e\t%.2e\t%.2e\t" ,dels,deldels,sigma_0,sud,logrr);
	//printf("%.2e\t%.2e\t%.2e\t%.2e\n" ,R,x, Q2,r);
	//getchar();
	val=(2.0-logrr)*dels+r*logrr*(dels*dels-deldels) ;
	//printf("%.2e\t" ,val);
	val*=exp(-sud)*BASE_SIGMA(r,x,Q2,SIGPAR);
	//printf("%.2e\n" ,val);
	
	return val;
}




double integral_term(double r, double x, double Q2,const  double * sigmapar,const  double* sudpar){
	//clock_t time=clock();
	VAR[0]=r;
	VAR[1]=x;
	VAR[2]=Q2;
	SIGPAR=sigmapar;
	SUDPAR=sudpar;

	double result=0;
	double rmin=1.0e-6;
	
	
	//double mu2;//=rmu2(r, SUDPAR);
	//compute_mu2(r,SUDPAR, &mu2, 1);
	
	//printf("1/rmin^2 =%.1e 1/r^2=%.1e \n" ,invrmin2,1.0/(r*r));
	//if( mu2>Q2 ){ 
	//	return(0.0);
	//}
#if SUDAKOV<=1
	rmin= pow( rmin2(Q2, SUDPAR ) ,0.5);
#endif
	//printf("rmin =%.1e r=%.1e \n" ,rmin,r);
	
	////////////////////////////////////////////
	//int N=96;
	//result=dgquad_(&integrand,&rmin,VAR,&N);
	///////////////////////////////////////////
	//double N=DGAUSS_PREC*0.01;
	//result=dgauss_(&integrand,&rmin,VAR,&N);
	////////////////////////////////////
	int seg=1;
	double NRel=DGAUSS_PREC;
	double NAbs=1.0e-10;
	double error=0;
	dadapt_(&integrand,&rmin,VAR,&seg ,&NRel, &NAbs, &result, &error)	;
	/////////////////////////////////////////
	
	return result;
	
}

double sigma_s(double r, double x, double Q2, const double * sigmapar, const double* sudpar){
	double sud=0.0;
	//double mu2=rmu2(r, sudpar);
	double mu2;
	compute_mu2(r,sudpar, &mu2, 1);
	//printf("%f\n",(&mu2)[0]);
	
	double val=BASE_SIGMA(r,x,Q2, sigmapar );
	//printf("%.3e\t%.3e\t%.3e\n",sudpar[0],pow(r,-2.0),pow(sudpar[1],-2.0));
#if SUDAKOV==0
	return val;
#endif
	//printf("1/rmin^2 =%.1e 1/r^2=%.1e \n" ,invrmin2,1.0/(r*r));
	//if( mu2>Q2 ){ 
	//	return(val);
	//}
	
#if SUDAKOV>=1
	
	sud+=sudakov_p(mu2,Q2,sudpar);
#endif
#if SUDAKOV==2
	sud+=sudakov_np(r,mu2,Q2,sudpar);
#endif
	//double sud=exp(-sudakov(r,Q2, sudpar ) );
	val*=exp(- sud);
	

	val+=integral_term(r,x,Q2,sigmapar,sudpar);
	//printf("%.2e \n",val);
	return val;
}





