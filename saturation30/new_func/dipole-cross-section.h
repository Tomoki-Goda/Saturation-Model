#define BGK 0

extern  double xgpdf(double, double);
extern double dgquad_(double (*)(double*), double*,double*,int*  );
////////////////////////////////////////////////////////////
////////////////// common functions ////////////////////////
////////////////////////////////////////////////////////////

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
//in simps.f
//extern "C" double simps_(double *, double *, double *,double *,double* ,double(*)(double*), double *  ,double *,double* ,double *);

/////////////////////////////////////////////////////////////
//////////////////////////// GBW ////////////////////////////
/////////////////////////////////////////////////////////////

double sigma_gbw(double r,double x,double Q2, double * par){
	double sigma_0 =par[0];
	double lambda	=par[1];
	double x_0	=par[2];
	
	//double xm=mod_x(x);

	return( sigma_0*(1-exp( - pow(r * Q0, 2) * pow(x_0/x, lambda)/4)) );	
}


/////////////////////////////////////////////////////////////
//////////////////////////// BGK ////////////////////////////
/////////////////////////////////////////////////////////////
#if BGK

//extern "C" double xgpdf(double, double);

double sigma_bgk(double r, double x, double Q2, double * par){
	double sigma_0		=par[0];
	double A_g		=par[1];
	double lambda_g		=par[2];
	double C		=par[3];
	double mu02		=par[4];
	
	//double xm=mod_x(x);	
	double mu2=C/(r*r)+mu02;
	double expo = (pow(r * PI,2) * alpha_s(mu2)* xgpdf(x,mu2))/ (3* sigma_0);
	
	return( sigma_0*(1-exp(-expo))) ;	
}
#else
double sigma_bgk(double r, double x, double Q2, double * par){
	return(0.0);
};
#endif

/////////////////////////////////////////////////////////////
//////////////////////////// GBS ////////////////////////////
/////////////////////////////////////////////////////////////
double sudakov(double r, double mu2,double* par) {
	//pertubative+non-perturbative sudakov
	double C= (*(par));
	double r_max=( *(par+1));
	double g1=(*(par+2));
	double g2=(*(par+3));

	double mub2=C*( 1.0/(pow(r,2)) + 1.0/pow(r_max,2) ) ;
	
	if (mu2 < LQCD2 || mub2 < LQCD2|| mu2 < mub2) {
		return(0.0);
	}
	
	double b0 = (11*CA-2*NF)/12;
	double val = CA/(2*b0*PI)*(log(mu2/LQCD2)*log(log(mu2/LQCD2)/log(mub2/LQCD2))-log(mu2/mub2));
#if SUDAKOV>=2	
	val+=(g1/2) * pow(r,2) + (g2/4) * log(1+pow(r/r_max,2) )*log( mu2 /pow(Q0 ,2) );
#endif
    return val;
}

#if SIMPS_GBS==1
double integrand_gbs(double r, double *par[2] ){
	//for the integral integrand has to be in the form 
	//func(double , double**) where integration is over the first argument. second are constants to be passed.
	//hence the following constant  parameters.
	if(fabs(r)<1.0e-15){
		return(0.0);
	}
	double R=( *(*par) ) ;
	double x=( *(*par+1) );
	//double xm=mod_x(x);
	double Q2=( *(*par+2) ) ;
	double sigma_0=(*(*(par+1) ));
	double lambda=( *(*(par+1)+1) );
	double x_0   =( *(*(par+1) +2) );
	double *sudpar=( *(par+1)+3 );//whatever parameter sudakov takes...

	double Qs2 =pow(Q0,2)*pow(x_0/x, lambda);
	double laplacian_sigma=sigma_0*r *log(R/r)*exp(-Qs2*pow(r,2) /4)*Qs2*(1-(Qs2*pow(r,2))/4);
	double val=laplacian_sigma;
#if SUDAKOV>=1
        val*=exp(-sudakov(r,Q2,sudpar)) ;
#endif
	//printf("integrand return for r=%f, %f. \n",r,val);       
	return(val); //sudakov(r,Q2,sudpar) *laplacian_sigma);
}


double sigma_gbs(double r, double x, double Q2, double * par){
	double param[3]={r,x,Q2};

	double *param_ptr[2]={param, par};
	double result=0;
	double error=0;
	simpson1dA(&integrand_gbs, param_ptr,0.0,r,200,&result,&error);
	//printf("\n\n");
	return(result);
}
#elif SIMPS_GBS==0
static double VAR[3];
static double PAR[5];

double integrand_gbs(double *r_ptr){
	//for the integral integrand has to be in the form 
	//func(double , double**) where integration is over the first argument. second are constants to be passed.
	//hence the following constant  parameters.
	double r=*r_ptr;
	if(fabs(r)<1.0e-15){
		return(0.0);
	}
	double R=( *(VAR)) ;
	double x=( *(VAR+1) );
	//double xm=mod_x(x);
	double Q2=( *(VAR+2) ) ;
	double sigma_0=(*(PAR ));
	double lambda=( *(PAR+1) );
	double x_0   =( *(PAR +2) );
	double *sudpar=( PAR+3 );//whatever parameter sudakov takes...

	double Qs2 =pow(Q0,2)*pow(x_0/x, lambda);
	double laplacian_sigma=sigma_0*r *log(R/r)*exp(-Qs2*pow(r,2) /4)*Qs2*(1-(Qs2*pow(r,2))/4);
#if SUDAKOV>=1
        val*=exp(-sudakov(r,Q2,sudpar)) ;
#endif
	//printf("integrand return for r=%f, %f. \n",r,val);       
	return(val); //sudakov(r,Q2,sudpar) *laplacian_sigma);
}

double sigma_gbs(double r, double x, double Q2, double * par){
	*(VAR)=r;
	*(VAR+1)=x;
	*(VAR+2)=Q2;
	
	for(unsigned i=0;i<5;i++){
		*(PAR+i)=*(par+i);
	}//not sure of efficiency..
		
	double result=0;
	double rmin=0.0;
	int N=96;
	
	result=dgquad_(&integrand_gbs,&rmin,VAR,&N);
	//printf("\n\n");
	return(result);
}
#endif

