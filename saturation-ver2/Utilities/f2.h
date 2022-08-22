



extern void simpson1dA(double(*)(double, double**),double**,double,double,int,double*,double*); 
extern double dgauss_(double (*)(const double*), double*,double*,double *  );
extern double dgquad_(double (*)(const double*),  const double*, const double*, const int*  );

///////////////////////////////////////Version 1////////////////////////////////////////////////////////
double f2_integrand(double r, double ** par){
	double xm;
	double x = *(*par);
	double Q2= *(*(par)+1);
	double *sigpar=*(par+1);
	double *sudpar=*(par+2);
	int fl=(int)(**( par+3)+0.1);

	//parameter(*(par+1),sigpar,sudpar);
	//printf("%f\t%f\t%f\n",par[1][0],par[1][1],par[1][2]);

	double value=0.0;

	//for(unsigned fl=0;fl<1/* NF-1*/;fl++){
//	for(unsigned fl=0;fl< 1;fl++){
		//printf("%d",fl);
		xm=mod_x(x,Q2,fl);
		value+=psisq_z_int(r, Q2, fl)* SIGMA(r,xm,Q2, sigpar,sudpar)/r ;
		
	//}
	//printf("x: %.2e\tx_mod %.2e\t Q2: %.2e\t%f\t%f\t%f\tresult :%f \n", x,mod_x(x,Q2,3),Q2,*(sigpar),*(sigpar+1),*(sigpar+2),value);
	return(value);
}
double f2(double Q2, double**pars){
	double res=0, err=0,val=0;
	*(*(pars)+1)=Q2;
	//printf("x== %.2e\t Q2== %.2e\n", pars[0][0],pars[0][1]);
	for(int i =0;i<(NF-1);i++){
		pars[3][0]=i;
		res=0;
		simpson1dA(&f2_integrand,pars,1.0e-5,30,250,&res,&err);
		val+=res;
	}
	//printf("%f\n",res);
	return val;
}

///////////////////////////////////////Version 2////////////////////////////////////////////////////////

static double Q2, X, *SUDPAR, *SIGPAR;
static int FL;

double f2_integrand_2(const double* R){
	double xm;

	double value=0;


	xm=mod_x(X,Q2,FL);
	double r=*R;
	value=psisq_z_int(r, Q2, FL)* SIGMA(r,xm,Q2,SIGPAR,SUDPAR)/(r) ;
		
	//printf("val=%.5e r=%.5e q2=%.5e  xm=%.5e\n",value,r ,Q2,xm );
	return(value);
}
double f2_2(double x,double q2, double *sigpar ,  double *sudpar){
	double res=0, err=0,val=0;
	SUDPAR=sudpar;
	SIGPAR=sigpar;
	//double prec=
	X=x;
	Q2=q2;
	int n=96;
	double rmin=R_MIN;
	double rmax=R_MAX;
	for(int i =0;i<(NF-1);i++){
		FL=i;
		res=0;
		//simpson1dA(&f2_integrand,pars,1.0e-5,30,250,&res,&err);
		res=dgquad_(&f2_integrand_2,&rmin,&rmax,&n);
		val+=res;
	}
	//printf("%f \tx=%.5e q2=%.5e fl=%d\n",val,X,Q2,FL);
	return val;
}










