#include<stdio.h>
#include<math.h>

#include"control.h"
#include"control-default.h"
#include"constants.h"

#include"./gluon-chebyshev.h"
#include"./sudakov.h"

extern double dgquad_(double (*)(double*), double*,double*,int*  );
extern double dgauss_(double (*)(double*), double*,double*,double *  );


////////////////////////////////////////////////////////////
////////////////// common functions ////////////////////////
////////////////////////////////////////////////////////////
int parameter(const double *par,double* sigpar,double* sudpar){
	sigpar[0]=par[0];
	sigpar[1]=par[1];
	sigpar[2]=par[2];
#if(MODEL==1||MODEL==3)
	sigpar[3]=par[3];
#if MU0==0
	sigpar[4]=par[4];
#else
	sigpar[4]=sqrt(fabs(sigpar[3]/par[4]));//rmax^2= C/mu02
#endif
#endif
////////////////////////////SUDPAR////////////////////////////////
#if (MODEL==22||MODEL==2)
	//for(int i=0;i<(N_PAR-3);i++){
	//	*(sudpar+i)=*(par+3+i);
	//}
	sudpar[0]=par[3];
#if MU0==1
	sudpar[1]=sqrt(fabs(sudpar[0]/par[4]));
#else
	sudpar[1]=par[4];
#endif
	
///////////////////////////////////////////////////////
#elif (MODEL==3)
	sudpar[0]=par[5] ;//C*C2
#if MU0==0 //if rmax is fit parameter
	sudpar[1]=par[4];//fabs is just in case;
#else //if mu02 is the fit parameter
	sudpar[1]=sqrt(fabs(sudpar[0]/par[4]));
#endif
#endif
//////////////////////////////////////////////////////
#if (SUDAKOV==2)
	sudpar[2]=par[6];
	sudpar[3]=par[7];
#endif
	return 0;
}


double alpha_s(double mu2 ){
	double b0= ((double)(33 -2*NF))/(12*PI);
	return( 1/(b0* log(mu2/LQCD2)));//LQCD2 lambda_QCD ^2
	
}
double mod_x(double x, double Q2, unsigned flavour) {
	double m_fsq;

	switch (flavour) {
	case 0:
		m_fsq = MASS_L2;
		break;
	case 1:
		m_fsq = MASS_S2;
		break;
	case 2:
		m_fsq = MASS_C2;
		break;
	case 3:
		m_fsq = MASS_B2;
		break;
	default:
		printf("mod_x::wrong input %c\n",flavour);
		m_fsq = MASS_L2;
	}
	return (x * (1.0 +( 4.0 * (m_fsq/Q2)) ));
}

/////////////////////////////////////////////////////////////
//////////////////////////// GBW ////////////////////////////
/////////////////////////////////////////////////////////////
double sigma_gbw(double r,double x,double q2, const double * par){
	double sigma_0 =par[0];
	double lambda	=par[1];
	double x_0	=par[2];
	
	if(x_0<0){//to avoid nan since migrad might give negative x0...
		return 0;
	}

	return( sigma_0*(1-exp( - pow(r * Q0, 2) * pow(x_0/x, lambda)/4)) );	
}



/////////////////////////////////////////////////////////////
//////////////////////////// BGK ////////////////////////////
/////////////////////////////////////////////////////////////
double sigma_bgk(double r, double x, double q2, const double * par){
	//clock_t tim=clock();
	double sigma_0		=par[0];
	double A_g		=par[1];
	double lambda_g	=par[2];
	//double C		=par[3];
	//double mu02		=par[4];
	//double rmax		=par[4];
	
	
	double mu2;
	compute_mu2(r, par+3 , &mu2, 1 );
	
	double expo = 0.389379*(pow( r* PI,2) * xg_chebyshev(x,mu2))/ (3* sigma_0); //prefactor, origin unknown...
	
	double val=sigma_0*(1-exp(-expo));
	
	return(val) ;	
}

//////////////////////////////////////Will be removed in the future/////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////             INTEGRATION            ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////  
#if (MODEL==2)

////////////////////////////////// CERN INTEGRATION ROUTINE VERSION ////////////////////////////////////
/// GLOBAL ///
static double VAR[3];
static const double *PAR;//its array but need only one because it needs only to point at par;

double integrand_gbs(double *r_ptr){
	double r=*r_ptr;
	if(fabs(r)<1.0e-15){
		return(0.0);
	}
	double R=( *(VAR)) ;
	double x=( *(VAR+1) );
	double Q2=( *(VAR+2) ) ;
	double laplacian_sigma=0;
	
	double sigma_0=(*(PAR ));
	double Qs2 =1;
#if (MODEL==2)
	double lambda=( *(PAR+1) );
	double x_0   =( *(PAR +2) );
	double *sudpar=( PAR+3 );//whatever parameter sudakov takes...
	 Qs2 =pow(Q0,2)*pow(x_0/x, lambda);
#endif

	laplacian_sigma=sigma_0*r *log(R/r)*exp(-Qs2*pow(r,2) /4)*Qs2*(1-(Qs2*pow(r,2))/4);


	double val=laplacian_sigma;
#if (SUDAKOV>=1)
        val*=exp(-sudakov(r,Q2,sudpar)) ;
#endif
	//printf("integrand return for r=%f, %f. \n",r,val);       
	return(val); //sudakov(r,Q2,sudpar) *laplacian_sigma);
}

double sigma_gbs(double r, double x, double Q2, const double * par){
	printf("discontinued\n");
	getchar();
	*(VAR)=r;
	*(VAR+1)=x;
	*(VAR+2)=Q2;
	 
	PAR=par;
		
	double result=0;
	double rmin=1.0e-5;
	double N=DGAUSS_PREC;
	//int N=96;
	//result=dgquad_(&integrand_gbs,&rmin,VAR,&N);
	result=dgauss_(&integrand_gbs,&rmin,VAR,&N);
	return(result);
}

#endif
