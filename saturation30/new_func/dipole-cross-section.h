//#define BGK 0

extern  double xgpdf(double, double);
extern double dgquad_(double (*)(double*), double*,double*,int*  );
extern double dgauss_(double (*)(double*), double*,double*,double *  );
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
double sigma_bgk(double r, double x, double Q2, double * par){
	//clock_t tim=clock();
	double sigma_0		=par[0];
	double A_g		=par[1];
	double lambda_g	=par[2];
	double C		=par[3];
	//double mu02		=par[4];
	double rmax		=par[4];
	
	
#if STAR==0
	double mu2=C*(1.0/(r*r)+1.0/(rmax*rmax)) ;
#elif STAR==1
	double mu02=C/(rmax*rmax);
	double mu2=mu02/(1-exp(-mu02 *pow(r,2)/C) );
#endif
	
	//double mu2=mu02/(1-exp(-pow(r/rmax,2)) );
	
	double expo = 0.389379*(pow( r* PI,2) * /*alpha_s(mu2)*/ xg_chebyshev(x,mu2))/ (3* sigma_0); //prefactor, origin unknown...
	
	double val=sigma_0*(1-exp(-expo));
	//tim-=clock();
	//printf("%.2e seconds\n",-((double)tim)/CLOCKS_PER_SEC);
	return(val) ;	
}


/////////////////////////////////////////////////////////////
//////////////////////////// GBS ////////////////////////////
/////////////////////////////////////////////////////////////
#if ((MODEL==2)||(MODEL==3)||(MODEL==22))
extern double rmu2( double ,double ,double );

double sudakov(double r, double mu2,double* par) {
	//pertubative+non-perturbative sudakov
	double C= (*(par));
	double r_max=( *(par+1));
	double g1=(*(par+2));
	double g2=(*(par+3));
	
	double mub2=rmu2(r ,C,r_max);
/*
#if STAR==0
	double mub2=C*(1.0/(r*r)+1.0/(r_max*r_max)) ;
#elif STAR==1
	double mu02=C/(r_max*r_max);
	double mub2=mu02/(1-exp(-pow(r/r_max,2) ) );
#endif
*/
	
	double val=0.0;

#if SUDAKOV>=2	
	val+=(g1/2) * pow(r,2) + (g2/4) * log(1+pow(r/r_max,2) )*log( mu2 /pow(Q0 ,2) );	
#endif
	
#if THETA_OFF==0
	if (mu2 < mub2) {//ensures that lower limit of integral is smaller than upper limit...
		return(val);
	}
#endif
		
//#if SUDAKOV>=2	
//	val+=(g1/2) * pow(r,2) + (g2/4) * log(1+pow(r/r_max,2) )*log( mu2 /pow(Q0 ,2) );	
//#endif
	
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
#endif
//////////////////////////////////////Will be removed in the future/////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////             INTEGRATION            ///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////  
#if MODEL==2
////////////////////////////////// CERN INTEGRATION ROUTINE VERSION ////////////////////////////////////
/// GLOBAL ///
static double VAR[3];
static double *PAR;//its array but need only one because it needs only to point at par;

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
/*	
#elif MODEL==3
	//double A_g=( *(PAR+1) );
	//double lambda_g  =( *(PAR +2) );
	double C =(*( PAR+3 ));
	double rmax =(*( PAR+4 ));
	double mu02=C/pow(rmax,2);
	double *sudpar=( PAR+3 );//whatever parameter sudakov takes...
	//double mu2=mu02/(1-exp(-pow(r/rmax,2)) );
	double mu2=C*(pow(r,-2)+pow(rmax,-2));
	//printf("s0 %.2e, A_g %.2e, l_g %.2e, C %.2e, rm %.2e ", *(PAR),*(PAR+1) ,*(PAR+2),*(PAR+3),*(PAR+4) );
	//printf("%.2e %.2e %.2e \n", R,x,Q2);
	
	double axg = xg_chebyshev(x,mu2);//\alpha_s(mu)x g(x,mu)...in chebyshev approx
	 Qs2 = 0.389379*4*PI*PI*axg/(3*sigma_0);
#endif	*/
	laplacian_sigma=sigma_0*r *log(R/r)*exp(-Qs2*pow(r,2) /4)*Qs2*(1-(Qs2*pow(r,2))/4);


	double val=laplacian_sigma;
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
	 /*for(unsigned i=0;i<N_PAR;i++){
                printf("%.3e ",*(par+i));
        }
        printf("\n");*/

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
