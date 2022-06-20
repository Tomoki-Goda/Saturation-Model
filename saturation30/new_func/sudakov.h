///////////////////////////////////////////////////////////////////////////////////////////////
extern double dgauss_(double(* )(double*),double*,double*,double*);

 

/////////////////////////////////////////////////////////////////////////////////////////////////

double dsdr(double r,double Q2,double C,double rmax){
#if SUDAKOV==0
	return 0.0;
#endif
	double mu2=C*(pow(r,-2)+pow(rmax,-2));
	if(mu2>Q2){
		return(0.0);
	}
	double jac= -2*C/pow(r,3);
	double val= - jac* alpha_s(mu2)*log(Q2/mu2)/mu2;
	return(CA/(2*PI) * val);
}

double ddsdrdr(double r,double Q2,double C,double rmax){
#if SUDAKOV==0
	return 0.0;
#endif
	double mu2=C*(pow(r,-2)+pow(rmax,-2));
	//jacobians (d mu2/d r )^2 and (d/dr)^2 mu2
	double jac1= 6*C/pow(r,4);
	double jac2=pow(2* C/pow(r,3),2);
	
	double b0 = (11*CA-2*NF)/12;
	if(mu2>Q2){
		return(0.0);
	}
	double val;

	
	double logQmu=log(Q2/mu2);
	double al=alpha_s(mu2);
	
	//derivative s wrt mu2 with jacobians
	double val1	= - jac1 * al*logQmu/mu2;	
	double val2	=   jac2 * al/pow(mu2,2) * ( 1 + logQmu*(1+b0*al));
	
	val=val1+val2;
	return(CA/(2*PI) * val); 
}

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
	double mu2=C* (pow(r,-2)+pow(rmax,-2));
	if(Q2<mu2){
		return 0.0;
	}
	
	double dels=dsdr(r,Q2,C,rmax);
	//double dels2=ddsdrdr(r,Q2,C,rmax);
	//val=(2-log(R/r) )*dels+r*log(R/r)*(pow(dels,2)-dels2);
	double Qs2=pow((*(PAR+2))/x ,*(PAR+1)  );
	val=dels* (BASE_SIGMA(r,x,Q2, PAR )*(1+ 2*log(R/r)-Qs2/4)/*+log(R/r) *exp(-r*r*Qs2/4  ) *(Qs2*r*r/2)  */);
	val*=exp(-sudakov(r,Q2, PAR+3 ) ) ;

	//val*=BASE_SIGMA(r,x,Q2, PAR );
	
	return val;	
}




double integral_term(double r, double x, double Q2,double* par){
	VAR[0]=r;
	VAR[1]=x;
	VAR[2]=Q2;
	PAR=par;
	
	double result=0;
	double rmin=1.0e-5;
	
	//double mu2=(*(par+3))*( pow(r,-2) +pow(*(par+4),-2) );
	rmin=pow(Q2/(*(par+3)) -pow(*(par+4),-2),-0.5);
	if( rmin>r ){
		printf(" no int term r=%.2e mu2=%e Q2=%.2e ",r ,rmin,Q2);	
	}
	
	double N=DGAUSS_PREC;
	//int N=96;
	//result=dgquad_(&integrand_gbs,&rmin,VAR,&N);
	result=dgauss_(&integrand,&rmin,VAR,&N);	
	
	return result;
	
}

double sigma_s(double r, double x, double Q2, double * par){
	double val=BASE_SIGMA(r,x,Q2, par );
	printf(" val=%.3e\t",val);
	double sud=exp(-sudakov(r,Q2, par+3 ) );
	
	val*=sud;
	double intterm=integral_term(r,x,Q2,par);
	printf("sud=%.3e\tintterm=%.3e\n",sud,intterm );
	
	val+=intterm;
	return val;
}
































