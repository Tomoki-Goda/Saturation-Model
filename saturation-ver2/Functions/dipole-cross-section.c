#include<stdio.h>
#include<math.h>

#include"control.h"
#include"control-default.h"
#include"constants.h"

#include"./gluon-chebyshev.h"
#include"./sudakov.h"
//#include"./laplacian.h"

extern double dgquad_(double (*)(double*), double*,double*,int*  );
extern double dgauss_(double (*)(double*), double*,double*,double *  );
extern double dclenshaw(double(*)(double*),double,double,double);

////////////////////////////////////////////////////////////
////////////////// common functions ////////////////////////
////////////////////////////////////////////////////////////
int parameter(const double *par,double* sigpar,double* sudpar){
///////////////////
//Sigpar are as we all know it, parameters for dipole sigma.
//sudpar are {C , r_max, g1, g2} but parameters may be given in terms of mu02 (as in BGK), and C and r_max may be that of BGK.
//That's why this is so messy...
////////////////////
	//for(int i=0;i<N_PAR;i++){
	//	printf("%.2e ",par[i]);
	//}
	//printf("\n");
#if (MODEL==0||MODEL==2||MODEL==22)
	sigpar[0]=par[0];
	sigpar[1]=par[1];
	sigpar[2]=par[2]*1.0e-4;
#else
	sigpar[0]=par[0];
	sigpar[1]=par[1];
	//sigpar[1]=pow(fabs(par[1])*1.0e-4,par[2]);
	sigpar[2]=par[2];
#endif

	
	
	//printf("SIGMA: %.2e %.2e %.2e ",sigpar[0],sigpar[1],sigpar[2]);
#if(MODEL==1||MODEL==3)
	sigpar[3]=par[3];
	#if MU0==0
		sigpar[4]=sigpar[3]/(par[4]*par[4]);
	#else
		sigpar[4]=par[4];//sqrt(fabs(sigpar[3]/par[4]));//rmax^2= C/mu02
	#endif
	//printf(" %.2e %.2e ",sigpar[3],sigpar[4]);
#endif
////////////////////////////SUDPAR////////////////////////////////
#if (MODEL==22||MODEL==2)
		sudpar[0]=par[3];
	#if MU0==0
		sudpar[1]=sudpar[0]/(par[4]*par[4]);
	#else
		sudpar[1]=par[4];
	#endif
	//printf("\tSUDAKOV: %.2e %.2e ",sudpar[0],sudpar[1]);	
#if (SUDAKOV==2)
	sudpar[2]=par[5];
	sudpar[3]=par[6];
	//printf("%.2e %.2e ",sudpar[2],sudpar[3]);
#endif

///////////////////////////////////////////////////////
#elif (MODEL==3)
	#if INDEPENDENT_C==1
		sudpar[0]=par[5] ;
	#else 
		sudpar[0]=par[3];
	#endif

	#if MU0==0 //if rmax is fit parameter
		#if INDEPENDENT_RMAX==1
			sudpar[1]=sudpar[0]/(par[6]*par[6]);
		#else
			sudpar[1]=sudpar[0]/(par[4]*par[4]);
		#endif

	#else //if mu02 is the fit parameter
		#if INDEPENDENT_RMAX==1
			sudpar[1]=par[6];
		#else
			sudpar[1]=par[4];//mu02 is shared
		#endif
	#endif
	//printf("%.2e %.2e ",sudpar[0],sudpar[1]);

//////////////////////////////////////////////////////
#if (SUDAKOV==2)
	sudpar[2]=par[7];
	sudpar[3]=par[8];
	//printf("%.2e %.2e ",sudpar[2],sudpar[3]);
#endif
#endif
	//printf("\n");
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
	//double A_g		=par[1];
	//double lambda_g		=par[2];
	//double C		=par[3];
	//double mu02		=par[4];
	//double rmax		=par[4];
	
	if(par[3]<0||par[4]<0){
		return 0;
	}
	double mu2;
	int signal= compute_mu2(r, par+3 , &mu2, 1 );
	if(signal!=0){
		printf("sigma_bgk:: C %.3e mu02 %.3e\n",par[3],par[4]);
		getchar();
		return 0;
	} 
	
	double expo = 0.389379*(pow( r* PI,2) * xg_chebyshev(x,mu2))/ (3* sigma_0); //prefactor, origin unknown...
	
	double val=sigma_0*(1-exp(-expo));
	
	return(val) ;	
}


////////////////////////////////////////////////////////////////////////////////////////
//derivatives
///////////////////////////////////////////////////////////////////////////////////////
struct parameters{
	double rmax;
	double x;
	double Q2;
	double *sigpar;
	double *sudpar;
} SIGPARAM;

void set_sigmapar(double x,double Q2 ,double*sigpar){
	SIGPARAM.x=x;
	SIGPARAM.Q2=Q2;
	SIGPARAM.sigpar=sigpar;
}

double sigma_for_lap(double *r){
	//double x=PARAM.x;
	double val=BASE_SIGMA(*r,SIGPARAM.x,1,SIGPARAM.sigpar);
	//double val=sigma_gbw(*r,SIGPARAM.x,1,SIGPARAM.sigpar);
	//printf("s=%.5e  r=%.5e\n",val,*r);
	return val;
}
double laplacian2(double (*func)(double*), double r,double step ){
	double arr[5];
	double x[5]={r-step,r-step/2,r,r+step/2,r+step};
	for(int i =0;i<5;i++){
		arr[i]=(*func)(x+i);
	}
	double d1=(-arr[4]+8*arr[3]-8*arr[1]+arr[0])/(6*step);
	double d2=(-arr[4]+16*arr[3]-30*arr[2]+16*arr[1]-arr[0])/(3*step*step);

	double val=d2+d1/r;//(arr[2]-2*arr[1]+arr[0])/(step*step)+(arr[2]-arr[0])/(r*step);
	return val;
}

double Qs2(double r , double x, double *par){
#if(MODEL==0||MODEL==2||MODEL==22)
	double qs2=pow(par[2]/x,par[1]);
#else
	double mu2;
	int signal= compute_mu2(r, par+3 , &mu2, 1 );
	double qs2 =4* 0.389379*(pow(  PI,2) * xg_chebyshev(x,mu2))/ (3* par[0]); //prefactor, origin unknown...
#endif
	return qs2;
}


double laplacian_sigma(double r,double x,double q2, double *par,double *sudpar){
	SIGPARAM.x=x;
	SIGPARAM.sigpar=par;
	double val;
	double bound=7.5e-3;
	if(r<bound){
		val=laplacian2(&sigma_for_lap,bound,bound/10);
		//val=laplacian2(&sigma_for_lap,r,r/5);
	}else if(r>3){
		double qs2=Qs2(r,x,par);

		val=par[0]*(1-r*r*qs2/4)*qs2*exp(-r*r*qs2/4);//(par[0]-BASE_SIGMA(r,x,1,par));
	}else{
		val=laplacian2(&sigma_for_lap,r,bound/10);
	}
#if SUDAKOV>=1
	double mu2;
	int signal=compute_mu2(r,sudpar, &mu2,1);//compute mu2
	if(q2>mu2){
		val*=exp_sud(r,mu2,q2);
	}
#endif
	return val;
}






//////////////////////////////////////Will be removed in the future/////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////             INTEGRATION            ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////  
double integrand_gbs(double *R){
	double r=*R;
	double val=r*log(SIGPARAM.rmax/r)*laplacian_sigma(r,SIGPARAM.x,SIGPARAM.Q2,SIGPARAM.sigpar,SIGPARAM.sudpar);
	return val;
}
double integrand_gbs2(double *R){
	double r=*R;
	double val=r*laplacian_sigma(r,SIGPARAM.x,SIGPARAM.Q2,SIGPARAM.sigpar,SIGPARAM.sudpar);
	return val;
}

double sigma_gbs(double r, double x, double Q2,double *sigpar,double* sudpar){
	set_sigmapar(x,Q2,sigpar);
	SIGPARAM.sudpar=sudpar;
	SIGPARAM.rmax=r;

	double rmin_2;
	int signal=rmin2(Q2, sudpar,&rmin_2 );

	double rmin=sqrt(rmin_2);
	double val=0;

	if(r-rmin<1.0e-15){
		val=BASE_SIGMA(r,x,Q2,sigpar);
		return(val);	
	}
	val=BASE_SIGMA(rmin,x,Q2,sigpar);
	double min=R_MIN;
	double eps=DGAUSS_PREC;
	val+=dclenshaw(&integrand_gbs,rmin,r,eps);
	val+=log(r/rmin)*dclenshaw(&integrand_gbs2,min,rmin,eps);
	printf("%.5e\n",val);
	return(val);
}


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
