///////////////////////////////////////////////////////////////////////////////////////////////
extern double dgauss_(double(* )(double*),double*,double*,double*); 
extern double dgquad_(double(* )(double*),double*,double*,int*);
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
	//clock_t time=clock();
	
	VAR[0]=r;
	VAR[1]=x;
	VAR[2]=Q2;
	PAR=par;

	double result=0;
	double rmin=1.0e-5;
	
	double mu2=(*(par+3))*( pow(r,-2) +pow(*(par+4),-2) );
	double invrmin2=Q2/(*(par+3)) - pow(*(par+4),-2);
	
	//printf("1/rmin^2 =%.1e 1/r^2=%.1e \n" ,invrmin2,1.0/(r*r));
	if( 1.0/(r*r)>invrmin2){ 
		return(0.0);
	}
	double N=DGAUSS_PREC;
	//int N=96;
	rmin=1.0/pow(invrmin2,0.5);
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
	//printf("first= %.3e sud=%.3e\t",val, sud);
	//double intterm=integral_term(r,x,Q2,par);
	//printf("intterm=%.3e ",intterm );
	//val+=intterm;
	val+=integral_term(r,x,Q2,par);
	//printf("\ttotal =%.3e\n",val);
	//for(unsigned i=0;i<N_PAR;i++){
	//	printf("%.2e ",par[i]);	
	//}
	//printf( "r=%.2e , x=%.2e Q2= %.2e \n",r,x,Q2);
	return val;
}





