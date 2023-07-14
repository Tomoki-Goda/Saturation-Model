#include"r-formula.hh"
//double Sigma::Qs2(const double x,const double r){
double Sigma_GBW::Qs2(const double x,const double r)const{
	if(x>1){
		printf("Sigma:: x is too large %.3e\n",x);
		getchar();
	}
	const double qs2=pow(x_0/x,lambda);//*pow(1-x,5.6); 
	if(qs2<1.0e-50){
		return 1.0e-50;
	}
	return qs2;		
}
double Sigma_BGK::Qs2(const double x,const double r)const{
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
	const double mu2=(rrmax>1.0e-6)?(mu02/((1.0-exp(-mu02*rrmax) ))):(1/rrmax+mu02/2+pow(mu02,2)*rrmax/12) ;
	//const double mu2=mu02/(1.0-exp(-mu02*rrmax)) ;
	const double al=alpha(mu2);
	//const double al=0.2;
	const double qs2=4*PI*PI*al*xgpdf(x,mu2,A_g,lambda_g)/(3*sigma_0); 
	if(qs2<1.0e-50){
		return 1.0e-50;
	}
	return qs2;		
}

inline double Sigma::freeze_x(const double x){
	double x2;
#if FREEZE_QS2==1
		x2=xmax*x/(xmax*(1-x)+x);
#elif  FREEZE_QS2==2
		x2=((x>xmax)?xmax:x);
#else
		x2=x;
#endif
/*
	switch(FREEZE_QS2){
	case 1:
		x2=0.5*x/(0.5*(1-x)+x);
		break;
	case 2:
		x2=((x>0.5)?0.5:x);
		break;
	default:
		x2=x;
	}
*/
	return x2;
}
//void Sigma::set_x(const double &x){
void Sigma_BGK::set_x(const double &x){
#if SIGMA_APPROX<0
	x2=freeze_x(x);
	xgpdf.set_x(x2);
#endif
}

//void Sigma::init(const double * const &sigpar){
void Sigma_GBW::init(const double * const &sigpar){
	this->par=sigpar;
	int i=0;
	this->sigma_0=sigpar[i++];
	this->lambda=sigpar[i++];
	this->x_0=sigpar[i++]; 
	#if MU02==0
	i++;//mu102=sigpar[i++];
	#endif
}
void Sigma_BGK::init(const double * const &sigpar){
	this->par=sigpar;
	int i=0;
	this->sigma_0=sigpar[i++];
	this->A_g=sigpar[i++];
	this->lambda_g=sigpar[i++];
	#if MU02==0
	i++;//mu102=sigpar[i++];
	#endif
	this->C=sigpar[i++];
	this->mu02=sigpar[i++];
	//printf("C=%.3e, mu02=%.3e, R_MIN=%.3e -> %.3e \n",C,mu02,R_MIN,C/pow(R_MIN,2)+mu02);
	#if SIGMA_APPROX<0
	this->xgpdf.init(mu02*0.5 , 2*(C/pow(R_MIN,2)+mu02),A_g,lambda_g);
	#endif

}

//inline double operator()(const double r)const {
//	return ((*this)(r,this->x));
//}
//double Sigma::operator()(const double x, const double r) {//,const double Q2,const double*sigpar)const {
double Sigma::operator()(const double x, const double r) {//,const double Q2,const double*sigpar)const {
	///Beware this transformation is also required in set_x()!!!!
	double x1=freeze_x(x);
	
	double qs2=Qs2(x1,r);
	#if ADJOINT==1
	qs2*=9.0/4.0;
	#endif
//	#if LAPLACIAN==0
	double val,ex;
	ex=pow(r,2)*qs2/4;
	if(ex<5.0e-3){
		//val*=(1-val/2+pow(val,2)/6-pow(val,3)/24);
		val=ex*(1-ex*(1-ex*(1-ex*(1-ex/5)/4)/3)/2);
		if(val>2.0e-3){
			if(fabs((1-exp(-ex))-val)>1.0e-15){
				printf("sigma discontinuous, %.2e %.2e \t %.3e\n",1-exp(-ex),val,val-(1-exp(-ex)));	
			}
		}
//		//printf("%.4e %.4e %.4e\n",r,(1-exp(-val))/pow(r,2),val/pow(r,2));
	}else{
		val=(1-exp(-ex));
	}
	val*=sigma_0;

	if(val<=0.0){
		printf("dp negative\n");
		val=0;
	}
	if(!std::isfinite(val)){
		return(0);
	}
	return val;
}
double Sigma:: S(const double x, const double r) {
	// S=sigma_0(1-sigma/sigma_0)
	//,const double Q2,const double*sigpar)const {
	///Beware this transformation is also required in set_x()!!!!
	double x1=freeze_x(x);
	
	double qs2=Qs2(x1,r);
	#if ADJOINT==1
	qs2*=9.0/4.0;
	#endif
	double val=sigma_0*exp(-pow(r,2)*qs2/4);

	if(!std::isfinite(val)){
		return(0);
	}
	return val;
}

