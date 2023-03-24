////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//   GBW / BGK dipoles
//   
//   Usage::
//   Sigma sigma // for BGK use Sigma<  (Cillinear_Gluon or Interpolate_Collinear_Gluon) >
// !! init for Collinear Gluon has to be fixed...
//   sigma.init(sigpar) //double* sigpar;
//   sigma(x,r) // double x,r 
//
// SIGMA 
////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if MODEL==1
template <typename ColG> class Sigma{
#else 
class Sigma{
#endif
// par and indivisual parameters are redundant. 
		//Collinear_Gluon xgpdf;
#if MODEL==1
		ColG xgpdf;
		double x2;
#endif
		double sigma_0,lambda, x_0, A_g,lambda_g,C,mu02,mu102,thresh_power;
		const double *par;
		inline double alpha(double mu2 ){
			static double b0= ((double)(33 -2*NF))/(12*PI);
			return( 1/(b0* log(mu2/LQCD2)));//LQCD2 lambda_QCD ^2
		}
		double Qs2(const double x,const double r){
			if(x>1){
				printf("x is too large %.3e\n",x);
				getchar();
			}
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
			const double rrmax=pow(r,2)/C;
			const double mu2=(r>1.0e-5)?(mu02/((1.0-exp(-mu02*rrmax) ))):(1/rrmax) ;
#if VARIANT<=2      
			const double al=alpha(mu2);
#elif VARIANT==3
			const double al=0.2;
#endif		
	
#if VARIANT==1
			const double qs2=4*PI*PI*al*A_g*pow(x,-lambda_g)*pow(1-x,5.6)/(3*sigma_0); 
#elif VARIANT==2
			const double qs2=4*PI*PI*al*A_g*pow(x,-lambda_g)/(3*sigma_0); 
#else
	#if SIGMA_APPROX>=0
		 	const double qs2=4*PI*PI*al*xgpdf(x,mu2,A_g,lambda_g)/(3*sigma_0); 
	#elif SIGMA_APPROX<0
			const double qs2=4*PI*PI*al*xgpdf(x,mu2)/(3*sigma_0); 
//		const double qs2=4*PI*PI*al*xgpdf(mu2)/(3*sigma_0); 
	#endif
#endif

#endif	//MODEL	
			return qs2;		
		}
	public:
		Sigma& operator=(const Sigma& rhs){
			init(rhs.par);
			return *this;
		} 
		explicit Sigma(void){ 
		}
		~Sigma(){
		}
		void set_x(const double &x){

#if FREEZE_QS2==1
			x2=0.5*x/(0.5*(1-x)+x);
#elif FREEZE_QS2==0
			x2=x;
#elif FREEZE_QS2==2
			x2=((x>0.5)?0.5:x);
#endif
#if SIGMA_APPROX<0
			xgpdf.set_x(x2);
#endif
		}
		void init(const double * const &sigpar){
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
			//printf("C=%.3e, mu02=%.3e, R_MIN=%.3e -> %.3e \n",C,mu02,R_MIN,C/pow(R_MIN,2)+mu02);
	#if SIGMA_APPROX<0
			xgpdf.init(5.0e-9,1,mu02*0.5,2*(C/pow(R_MIN,2)+mu02),A_g,lambda_g);
	#endif
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
		double operator()(const double x, const double r) {//,const double Q2,const double*sigpar)const {
#if FREEZE_QS2==1 ////Beware this transformation is also required in set_x()!!!!
			double x1=0.5*x/(0.5*(1-x)+x);
#elif FREEZE_QS2==0
			double x1=x;
#elif FREEZE_QS2==2
			double x1=((x>0.5)?0.5:x);
#endif///////////////////////////////////////////////////////////////////////////
			if(x1!=x2){
				printf("Sigma:: Error: x does not match. input=%.3e internal x= %.3e diff = %.3e\n",x1,x2, x1-x2);
			}
			double qs2=Qs2(x1,r);
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
			val*=pow(1-x1,thresh_power);
#endif
			if(!std::isfinite(val)){
				return(0);
			}
			return val;
		}

};
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////



