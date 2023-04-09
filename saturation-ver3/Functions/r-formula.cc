

#include"r-formula.hh"
//double Sigma::Qs2(const double x,const double r){
double Sigma_GBW::Qs2(const double x,const double r){
	if(x>1){
		printf("Sigma:: x is too large %.3e\n",x);
		getchar();
	}
	const double qs2=pow(x_0/x,lambda);//*pow(1-x,5.6); 
	if(qs2<1.0e-25){
		return 1.0e-25;
	}
	return qs2;		
}
double Sigma_BGK::Qs2(const double x,const double r){
	if(x>1){
		printf("Sigma:: x is too large %.3e\n",x);
		getchar();
	}
	#if SIGMA_APPROX<0
	if(x!=x2){
		printf("Sigma:: Error: x does not match. input=%.3e internal x= %.3e diff = %.3e\n",x,x2, x-x2);
	}
	#endif
	const double rrmax=pow(r,2)/C;
	const double mu2=(r>1.0e-5)?(mu02/((1.0-exp(-mu02*rrmax) ))):(1/rrmax) ;
	const double al=alpha(mu2);
	#if SIGMA_APPROX>=0
	const double qs2=4*PI*PI*al*xgpdf(x,mu2,A_g,lambda_g)/(3*sigma_0); 
	#elif SIGMA_APPROX<0
	const double qs2=4*PI*PI*al*xgpdf(x,mu2)/(3*sigma_0); 
	#endif
	if(qs2<1.0e-25){
		return 1.0e-25;
	}
	//printf("Qs2=%.2e\n",qs2);
	return qs2;		
}

//void Sigma::set_x(const double &x){
void Sigma_BGK::set_x(const double &x){
	#if SIGMA_APPROX<0
	#if FREEZE_QS2==1
	x2=0.5*x/(0.5*(1-x)+x);
	#elif FREEZE_QS2==0
	x2=x;
	#elif FREEZE_QS2==2
	x2=((x>0.5)?0.5:x);
	#endif
	xgpdf.set_x(x2);
	#endif
}

//void Sigma::init(const double * const &sigpar){
void Sigma_GBW::init(const double * const &sigpar){
	par=sigpar;
	int i=0;
	sigma_0=sigpar[i++];
	lambda=sigpar[i++];
	x_0=sigpar[i++]; 
	#if MU02==0
	i++;//mu102=sigpar[i++];
	#endif
	#if THRESHOLD==-1
	thresh_power=sigpar[i++];
	#else
	thresh_power=THRESHOLD;
	#endif
	//sigpar=par;
}
void Sigma_BGK::init(const double * const &sigpar){
	par=sigpar;
	int i=0;
	sigma_0=sigpar[i++];
	A_g=sigpar[i++];
	lambda_g=sigpar[i++];
	#if MU02==0
	i++;//mu102=sigpar[i++];
	#endif
	C=sigpar[i++];
	mu02=sigpar[i++];
	//printf("C=%.3e, mu02=%.3e, R_MIN=%.3e -> %.3e \n",C,mu02,R_MIN,C/pow(R_MIN,2)+mu02);
	#if SIGMA_APPROX<0
	xgpdf.init(mu02*0.5 , 2*(C/pow(R_MIN,2)+mu02),A_g,lambda_g);
	#endif

	#if THRESHOLD==-1
	thresh_power=sigpar[i++];
	#else
	thresh_power=THRESHOLD;
	#endif
	//sigpar=par;
}

//inline double operator()(const double r)const {
//	return ((*this)(r,this->x));
//}
//double Sigma::operator()(const double x, const double r) {//,const double Q2,const double*sigpar)const {
double Sigma::operator()(const double x, const double r) {//,const double Q2,const double*sigpar)const {
	#if FREEZE_QS2==1 ////Beware this transformation is also required in set_x()!!!!
	double x1=0.5*x/(0.5*(1-x)+x);
	#elif FREEZE_QS2==0
	double x1=x;
	#elif FREEZE_QS2==2
	double x1=((x>0.5)?0.5:x);
	#endif///////////////////////////////////////////////////////////////////////////
	
	double qs2=Qs2(x1,r);
	#if ADJOINT==1
	qs2*=9.0/4.0;
	#endif

	#if IBP==2
	double val=-sigma_0*exp(-pow(r,2)*qs2/4);
	#else
//	#if LAPLACIAN==0
	double val;
	val=pow(r,2)*qs2/4;
	if(val<1.0e-4){
		val*=(1-val/2+pow(val,2)/6-pow(val,3)/24);
		//printf("%.4e %.4e %.4e\n",r,(1-exp(-val))/pow(r,2),val/pow(r,2));
	}else{
		val=(1-exp(-val));
	}
	val*=sigma_0;

	if(val<=0.0){
		printf("dp negative\n");
		val=0;
	}
//	#elif LAPLACIAN==1
//	double val=pow(r,2)*qs2/4;
//	val=sigma_0*qs2*(1-val)*exp(-val);
//	#endif//LAPLACIAN
	#endif//IBP==2
	#if THRESHOLD!=0 
	//double thresh_power=THRESHOLD;
	val*=pow(1-x1,thresh_power);
	#endif
	if(!std::isfinite(val)){
		return(0);
	}
	return val;
}
