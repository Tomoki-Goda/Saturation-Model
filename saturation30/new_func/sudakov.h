///////////////////////////////////////////////////////////////////////////////////////////////
extern double dgauss_(double(* )(double*),double*,double*,double*); 

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////Jacobian/////////////////////////////////////////
double dsdr(double r,double Q2,double C,double rmax){
#if SUDAKOV==0
	return 0.0;
#endif
	double mu2=C*(pow(r,-2.0)+pow(rmax,-2.0));
	if(mu2>Q2){
		return(0.0);
	}
	if(mu2<LQCD2){
		printf( "drds:: mu2 out of range");
	}

	double jac= -2.0*C/pow(r,3.0);
	double b0= (11.0*CA-2.0*NF)/12.0;
	double alpha =1.0/(b0* log( mu2/LQCD2) );
	double val= - jac* alpha*log(Q2/mu2)/mu2;
	return(((double)CA)/(2.0*PI) * val);
}

double ddsdrdr(double r,double Q2,double C,double rmax){
#if SUDAKOV==0
	return 0.0;
#endif
	double mu2=C*(pow(r,-2.0)+pow(rmax,-2.0));
	//jacobians (d mu2/d r )^2 and (d/dr)^2 mu2
	double jac1= 6.0*C/pow(r,4.0);
	double jac2=pow(2.0* C/pow(r,3.0),2.0);
	
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
	double val1	= - jac1 * alpha*logQmu/mu2;	
	double val2	=   jac2 * alpha/pow(mu2,2) * ( 1 + logQmu*(1+b0*alpha));
	
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
	double mu2=C* (pow(r,2)+pow(rmax,2))/(r*r*rmax*rmax);
	
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
	VAR[0]=r;
	VAR[1]=x;
	VAR[2]=Q2;
	PAR=par;

	double result=0;
	double rmin=1.0e-5;
	
	double mu2=(*(par+3))*( pow(r,-2) +pow(*(par+4),-2) );
	rmin=1.0/pow(Q2/(*(par+3)) -pow(*(par+4),-2),0.5);
	double N=DGAUSS_PREC;
	//int N=96;
	//result=dgquad_(&integrand_gbs,&rmin,VAR,&N);
	result=dgauss_(&integrand,&rmin,VAR,&N);	
	
	return result;
	
}

double sigma_s(double r, double x, double Q2, double * par){
	double sud=exp(-sudakov(r,Q2, par+3 ) );
	
	double val=sud*BASE_SIGMA(r,x,Q2, par );
	
	double intterm=integral_term(r,x,Q2,par);
	//printf("first= %.3e sud=%.3e\tintterm=%.3e\n",val, sud,intterm );
	val+=intterm;
	printf("=%.3e\n",val);
	return val;
}





