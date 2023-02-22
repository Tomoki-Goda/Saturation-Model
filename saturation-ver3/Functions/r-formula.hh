
//extern double change_var(double & var,double &  jac,const double min, const double max,const double c);
//extern "C" double xgpdf_(const double& x, const double& QQ,const double& A_g, const double& lambda_g );
//extern "C" double xgpdf_(const double* x, const double* QQ,const double* A_g, const double* lambda_g );


////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//   GBW / BGK dipoles
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Sigma{
		Collinear_Gluon xgpdf;
		
		inline double alpha(double mu2 )const{
			static double b0= ((double)(33 -2*NF))/(12*PI);
			return( 1/(b0* log(mu2/LQCD2)));//LQCD2 lambda_QCD ^2
		}
		//inline double alpha(const double  mu2)const{
		//	return 4.0/(9.0 *log( ((mu2>2*LQCD2)?(mu2):(2.0*LQCD2))/LQCD2));
		//}
		const double* sigpar;
		double x=0;
	public:
		inline void set_kinem(const double x){
			this->x=x;
		}
		inline void init(const double *par){
			sigpar=par;
		}
		explicit Sigma(void){ 
		}
		~Sigma(){
		}
		inline double operator()(const double r,const double x){//,const double Q2,const double*sigpar)const {
		 	this->x=x;
		 	return((*this)(r));
		 }

		double operator()(const double r)const {
			const double sigma_0=sigpar[0];
#if MODEL==0//GBW
			const double lambda=sigpar[1];
			const double x_0=sigpar[2];
			if(x_0<1.0e-5||x_0>1.0e-3){
				return 0;
			}
			if(lambda<0.05||lambda>0.95){
				return 0;
			}
			const double qs2=pow(x_0/x,lambda);//*pow(1-x,5.6); 
#elif MODEL==1//BGK
			//set_xg_parameter(sigpar[1],sigpar[2]);
			const double C=sigpar[3];
			const double mu02=sigpar[4];
			if(sigpar[1]<=0||C<=0||mu02<1){
				return 0;
			}
			//double mu2=C/(r*r)+mu02;
			const double exprrmax=exp(-r*r*(mu02/C));
	
			const double mu2=mu02/((1.0-exprrmax ));
			if(!std::isfinite(mu2)){
				printf("r= %.3e C=%.3e mu02=%.3e ,mu2=%.3e \n",r,C,mu02,mu2 );
				return 0;
				//getchar();

			}
			if(mu2<1||!std::isfinite(mu2)){
				return(0);
			}
			/*
#if FREEZE_QS2==1
			const double qs2=4*PI*PI*alpha(mu2)*xg(((x>0.5)?0.5:x),mu2,sigpar[1],sigpar[2])/(3*sigma_0); 
#else
			const double qs2=4*PI*PI*alpha(mu2)*xg(x,mu2,sigpar[1],sigpar[2])/(3*sigma_0); 
#endif*/
#if FREEZE_QS2==1
			const double qs2=4*PI*PI*alpha(mu2)*xgpdf(((x>0.5)?0.5:x),mu2,sigpar[1],sigpar[2])/(3*sigma_0); 
#else
			//const double qs2=4*PI*PI*alpha(mu2)*xgpdf_(&x,&mu2,sigpar+1,sigpar+2)/(3*sigma_0);
			const double qs2=4*PI*PI*alpha(mu2)*xgpdf(x,mu2,sigpar[1],sigpar[2])/(3*sigma_0); 
#endif
			//double qs2=4*PI*PI*alpha(mu2)*xgpdf_(&x,&mu2,sigpar+1,sigpar+2)/(3*sigma_0); 
			//double qs2=4*PI*PI*alpha(mu2)*xgpdf(x,mu2)/(3*sigma_0); 
#endif	//MODEL	

//#if (LAPLACIAN==0||R_FORMULA==1)
#if IBP==2
			double val=-sigma_0*exp(-pow(r,2)*qs2/4);
#else
#if LAPLACIAN==0
			double val;
			val=pow(r,2)*qs2/4;
			if(val<1.0e-6){
				val*=sigma_0;
			}else{
				val=sigma_0*(1-exp( -val));
			}
#elif LAPLACIAN==1
			double val=sigma_0*qs2*(1-r*r*qs2/4)*exp(-r*r*qs2/4);
#endif//LAPLACIAN
#endif//IBP==2
#if THRESHOLD==1 
			double thresh_power=7;
			val*=pow(1-x,thresh_power);
#endif
			//printf("%.3e\n",val);
			//printf("val=%.3e: r= %.3e C=%.3e mu02=%.3e ,mu2=%.3e, Qs2=%.3e \n",val, r,C,mu02,mu2,qs2 );
			if(!std::isfinite(val)){
				return(0);
			}
			return val;
		}

};
