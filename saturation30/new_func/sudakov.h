///////////////////////////////////////////////////////////////////////////////////////////////
extern double dgauss_(double(* )(double*),double*,double*,double*); 
extern double dgquad_(double(* )(double*),double*,double*,int*);
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////Something for the parametrizations of mu2/////////////////////////////////
///////////////////////////////Jacobian dmu2/dr and ddmu2/drdr/////////////////////////////////////////

double rmu2( double r ,double C,double rmax){

	double mu2;
#if STAR==0
	mu2=C*(pow(r,-2.0)+pow(rmax,-2.0));
#elif STAR==1 //new parametrization
	double exprrmax=exp(-pow(r/rmax,2));
	mu2=C/(rmax*rmax*(1.0-exprrmax ));
#endif
	return mu2;
}
//////////////////////////
double rmu2_jac_first( double r ,double C,double rmax){
	double jac;
#if STAR==0
	jac= -2.0*C/pow(r,3.0);
#elif STAR==1 //new parametrization  
	double exprrmax=exp(-pow(r/rmax,2));
	jac=-2.0*C*(r/pow(rmax,4))*( exprrmax/pow(1.0-exprrmax,2));
#endif
	return jac;
}
//////////////////////////
double rmu2_jac_second( double r ,double C,double rmax){
	double jac;
#if STAR==0
	jac=6.0*C/pow(r,4.0);
#elif STAR==1 //new parametrization  
	double exprrmax=exp(-pow(r/rmax,2));
	jac =- 2.0*(C/pow(rmax,4))*(exprrmax/pow(1.0-exprrmax,2))*( 1.0- 2* pow(r/rmax,2)*((1.0+exprrmax)/(1.0-exprrmax)) );
#endif
	return jac;
}

/////////////////////////////////////////////
double rmin2(double Q2,double C,double rmax ){
	//if (Q2<(C/(rmax*rmax))){
	//	return 0.0;//Q2 too small;
	//}
	double rmin2;
#if STAR==0
	rmin2=-rmax*rmax*C/(C-Q2*rmax*rmax);
#elif STAR==1
	rmin2=-rmax*rmax*log( 1.0-C/(Q2*rmax*rmax) );
#endif
	return rmin2;
}
///////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////jacobian ////////////////////////////////////////////

double dsdr(double r,double Q2,double C,double rmax){
#if SUDAKOV==0
	return 0.0;
#endif
	
	double mu2=rmu2(r,C,rmax);
	double jacfirst=rmu2_jac_first(r,C,rmax);
	
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

double ddsdrdr(double r,double Q2,double C,double rmax){
#if SUDAKOV==0
	return 0.0;
#endif
	double mu2=rmu2(r,C,rmax);
	double jacfirst=rmu2_jac_first(r,C,rmax);
	double jacsecond=rmu2_jac_second(r,C,rmax);
	
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
static double *PAR;


double integrand( double * r_ptr){
	double val;
	double r=*r_ptr;
	double  R=VAR[0];
	double  x=VAR[1];
	double  Q2=VAR[2];
	
	double C=PAR[3];
	double rmax=PAR[4];
	double mu2=rmu2(r,C,rmax);
	
	if(Q2<mu2){
		return 0.0;
	}
	if(mu2<LQCD2){
		printf("integrand::mu2 out of range");
	}
		
	double dels=dsdr(r,Q2,C,rmax);
	double deldels=ddsdrdr(r,Q2,C,rmax);

	double sigma_0=*(PAR);
	double sud=sudakov(r,Q2, PAR+3 );
	double logrr=log(R/r);
	
	val=(2.0-logrr)*dels+r*logrr*(dels*dels-deldels) ;
	val*=exp(-sud)*BASE_SIGMA(r,x,Q2,PAR);
}




double integral_term(double r, double x, double Q2,double* par){
	//clock_t time=clock();
	
	VAR[0]=r;
	VAR[1]=x;
	VAR[2]=Q2;
	PAR=par;

	double result=0;
	double rmin=1.0e-5;
	double C=*(par+3);
	double rmax=*(par+4);
	double mu2=rmu2(r,C,rmax);
	
	//printf("1/rmin^2 =%.1e 1/r^2=%.1e \n" ,invrmin2,1.0/(r*r));
	if( mu2>Q2 ){ 
		return(0.0);
	}
	rmin= pow( rmin2(Q2,C,rmax ) ,0.5);
	//printf("rmin =%.1e r=%.1e \n" ,rmin,r);
	double N=DGAUSS_PREC;
	//int N=96;
	//rmin=1.0/pow(invrmin2,0.5);
	//printf("rmin= %.1e\n",rmin);
	//result=dgquad_(&integrand,&rmin,VAR,&N);
	result=dgauss_(&integrand,&rmin,VAR,&N);	
	
	//time-=clock();
	//printf("%.2e seconds\n", -((double)time)/CLOCKS_PER_SEC);
	
	return result;
	
}

double sigma_s(double r, double x, double Q2, double * par){
	double sud=exp(-sudakov(r,Q2, par+3 ) );
	
	double val=sud*BASE_SIGMA(r,x,Q2, par );
	val+=integral_term(r,x,Q2,par);
	
	return val;
}





