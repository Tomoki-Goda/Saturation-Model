////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//   GBW / BGK dipoles
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Sigma{
// par and indivisual parameters are redundant. 
		Collinear_Gluon xgpdf;
		double x=0;
		double sigma_0,lambda, x_0, A_g,lambda_g,C,mu02,mu102,thresh_power;
		const double *par;
		
		//double Qs2(const double x,const double r,const double (&sigpar)[])const{
		inline double alpha(double mu2 )const{
			static double b0= ((double)(33 -2*NF))/(12*PI);
			return( 1/(b0* log(mu2/LQCD2)));//LQCD2 lambda_QCD ^2
		}
		double Qs2(const double x,const double r)const{
#if MODEL==0//GBW
#if VARIANT==0
			const double qs2=pow(x_0/x,lambda);//*pow(1-x,5.6); 
#elif VARIANT==1
			const double qs2=pow(x_0/x,lambda)*pow(1-x,5.6); 
#elif VARIANT==2
			const double exprrmax=exp(-pow(r,2)*(1/0.25));
			const double mu2=1/((1.0-exprrmax ));
			const double qs2=alpha(mu2)*pow(x_0/x,lambda);//*pow(1-x,5.6); 
#elif VARIANT==3
			const double exprrmax=exp(-pow(r,2)*(1/0.25));
			const double mu2=1/((1.0-exprrmax ));
			const double qs2=alpha(mu2)*pow(x_0/x,lambda)*pow(1-x,5.6); 
#endif
#elif MODEL==1//BGK
			const double exprrmax=exp(-pow(r,2)*(mu02/C));
			const double mu2=mu02/((1.0-exprrmax ));
#if FREEZE_QS2==1
			double x1=0.5*x/(0.5*(1-x)+x);
#elif FREEZE_QS2==0
			double x1=x;
#elif FREEZE_QS2==2
			double x1=((x>0.5)?0.5:x);
#endif

#if VARIANT<=2      
			const double al=alpha(mu2);
#elif VARIANT==3
			const double al=0.2;
#endif			
#if VARIANT==1
			const double qs2=4*PI*PI*al*A_g*pow(x,-lambda_g)*pow(1-x1,5.6)/(3*sigma_0); 
#elif VARIANT==2
			const double qs2=4*PI*PI*al*A_g*pow(x,-lambda_g)/(3*sigma_0); 
#else
			const double qs2=4*PI*PI*al*xgpdf(x1,mu2,A_g,lambda_g)/(3*sigma_0); 
#endif

#endif	//MODEL	
			return qs2;		
		}
		
		
	public:
		Sigma& operator=(const Sigma& rhs){
			x=rhs.x;
			init(rhs.par);
			return *this;
		} 
		explicit Sigma(const Sigma& rhs){
			//xgpdf=rhs.xgpdf;
			x=rhs.x;
			init(rhs.par);
		}
		explicit Sigma(void){ 
		}
		~Sigma(){
		}
		inline void set_kinem(const double x){
			this->x=x;
		}
		inline void init(const double * const &sigpar){
			par=sigpar;
			int i=0;
#if MODEL==0//GBW
			sigma_0=sigpar[i++];
			lambda=sigpar[i++];
			x_0=sigpar[i++]; 
#elif MODEL==1//BGK
			sigma_0=sigpar[i++];
			A_g=sigpar[i++];
			lambda_g=sigpar[i++];
#endif
#if MU02==0
			i++;//mu102=sigpar[i++];
#endif
#if MODEL==1
			C=sigpar[i++];
			mu02=sigpar[i++];
#endif

#if THRESHOLD==-1
			thresh_power=sigpar[i++];
#else
			thresh_power=THRESHOLD;
#endif
			//sigpar=par;
		}
		
		inline double operator()(const double r,const double x){//,const double Q2,const double*sigpar)const {
		 	this->x=x;
		 	return((*this)(r));
		 }

		double operator()(const double r)const {
			//const double sigma_0=sigpar[0];
			double qs2=Qs2(x,r);
#if IBP==2
			double val=-sigma_0*exp(-pow(r,2)*qs2/4);
#else
#if LAPLACIAN==0
			double val;
			val=pow(r,2)*qs2/4;
			if(val<1.0e-7){
				val*=sigma_0;
			}else{
				val=sigma_0*(1-exp( -val));
			}
#elif LAPLACIAN==1
			double val=pow(r,2)*qs2/4;
			val=sigma_0*qs2*(1-val)*exp(-val);
#endif//LAPLACIAN
#endif//IBP==2
#if THRESHOLD!=0 
			//double thresh_power=THRESHOLD;
			val*=pow(1-x,thresh_power);
#endif
			if(!std::isfinite(val)){
				return(0);
			}
			return val;
		}

};
