



extern void simpson1dA(double(*)(double, double**),double**,double,double,int,double*,double*); 
extern double dgauss_(double (*)(const double*), double*,double*,double *  );
extern double dgquad_(double (*)(const double*),  const double*, const double*, const int*  );
extern double dadapt_(double(* )(const double*),double*,double*,int*,double*,double* ,double*,double*);

static int FIX_W=0;
int PLOT_FLAVOUR=0;

double mod_x_W(double x, double Q2,double W2,int fl){
	double mass2;

	switch(fl){
		case 0:
			mass2=MASS_L2;
			break;
		case 1:
			mass2=MASS_S2;
			break;
		case 2:
			mass2=MASS_C2;
			break;
		case 3:
			mass2=MASS_B2;
			break;
		default:
			printf("psisq_f::wrong input %d\n",fl);
			mass2=MASS_L2;
	}
	
	double xm=(Q2+4*mass2)/(Q2+W2);
	return  xm;
}

///////////////////////////////////////Version 1////////////////////////////////////////////////////////
double f2_integrand(double R, double ** par){
	double xm;
	double x = *(*par);
	double Q2= *(*(par)+1);
	double *sigpar=*(par+1);
	double *sudpar=*(par+2);
	int fl=(int)(**( par+3)+0.1);
	
	double r=R/(1-R);
	double jac=pow(1-R,-2);

	double value=0.0;

	if(FIX_W==1){
		double W2=pow(*(*(par)),2);
		xm=mod_x_W(x, Q2, W2,fl);
	}else{
		xm=mod_x(x,Q2,fl);
	}
	value+=jac*psisq_z_int(r, Q2, fl)* SIGMA(r,xm,Q2, sigpar,sudpar)/r ;
	return(value);
}
double f2(double Q2, double**pars){
	double res=0, err=0,val=0;
	*(*(pars)+1)=Q2;
	printf("fl=%d\n",PLOT_FLAVOUR);

	//printf("x== %.2e\t Q2== %.2e\n", pars[0][0],pars[0][1]);
	if(PLOT_FLAVOUR==0){
		for(int i =0;i<(NF-1);i++){
			pars[3][0]=i;
			res=0;
			simpson1dA(&f2_integrand,pars,1.0e-5,0.97,500,&res,&err);
			val+=res;
		}
	}else if(PLOT_FLAVOUR==2){
		pars[3][0]=2;
		res=0;
		simpson1dA(&f2_integrand,pars,1.0e-5,0.97,500,&res,&err);
		val+=res;
	}else if(PLOT_FLAVOUR==3){
		pars[3][0]=3;
		res=0;
		simpson1dA(&f2_integrand,pars,1.0e-5,0.97,500,&res,&err);
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
	double r=*R/(1-*R);
	double jac=pow(1-*R,-2);
	value=jac*psisq_z_int(r, Q2, FL)* SIGMA(r,xm,Q2,SIGPAR,SUDPAR)/(r) ;
		
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
	//
	double rmin=R_MIN;
	double rmax=R_MAX;
	rmin=rmin/(1+rmin);
	rmax=rmax/(1+rmax);
	//int seg=30;
	//double NRel=1.0e-1;
	//double NRel=DGAUSS_PREC;
	//double NAbs=0;
	//double error=0;
	
	double N=DGAUSS_PREC;
	int n=96;
	double step=(rmax-rmin)/7;
	double high,low;
	
	if(PLOT_FLAVOUR==0){
		for(int i =0;i<(NF-1);i++){
			FL=i;
			res=0;
			low=rmin;
			high=rmin;
			//simpson1dA(&f2_integrand,pars,1.0e-5,30,250,&res,&err);
			//dadapt_(&f2_integrand_2,&rmin,&rmax,&seg ,&NRel, &NAbs, &res, &error);
			//res+=dgauss_(&f2_integrand_2,&rmin,&rmax,&N);
			//res+=dgquad_(&f2_integrand_2,&rmin,&rmax,&n);
			for(int i=0;i<7;i++){
				high+=step;
				//res+=dgauss_(&f2_integrand_2,&low,&high,&N);
				res+=dgquad_(&f2_integrand_2,&low,&high,&n);
				//printf(" %f %f %d\n",rmin,rmax,n);
				low=high;
			}
			//printf("%.5e %.5e\n",res,error);
			val+=res;
		}
	}else {
		FL=PLOT_FLAVOUR;
		res=0;
		low=rmin;
		high=rmin;
		//simpson1dA(&f2_integrand,pars,1.0e-5,30,250,&res,&err);
		//dadapt_(&f2_integrand_2,&rmin,&rmax,&seg ,&NRel, &NAbs, &res, &error);
		//res+=dgauss_(&f2_integrand_2,&rmin,&rmax,&N);
		//res+=dgquad_(&f2_integrand_2,&rmin,&rmax,&n);
		for(int i=0;i<7;i++){
			high+=step;
			//res+=dgauss_(&f2_integrand_2,&low,&high,&N);
			res+=dgquad_(&f2_integrand_2,&low,&high,&n);
			//printf(" %f %f %d\n",rmin,rmax,n);
			low=high;
		}
		//printf("%.5e %.5e\n",res,error);
		val+=res;
	}
	//printf("%f \tx=%.5e q2=%.5e fl=%d\n",val,X,Q2,FL);
	return val;
}










