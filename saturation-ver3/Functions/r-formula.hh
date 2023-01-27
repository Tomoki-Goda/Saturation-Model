#include"./gluons.hh"
inline PREC modx(const PREC x, const PREC Q2, const  PREC mf2){
#if MODX==1
	return( (x*(1+4*mf2/Q2)));
#else 
	return( x);
#endif
}
extern double change_var(double & var,double &  jac,const double min, const double max,const double c);
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//   GBW
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Sigma{
		Collinear_Gluon xg;
		
	//	inline double alpha(double mu2 ){
	//		static double b0= ((double)(33 -2*NF))/(12*PI);
	//		return( 1/(b0* log(mu2/LQCD2)));//LQCD2 lambda_QCD ^2
	//	}
		inline double alpha(const PREC mu2)const{
			return 4.0/(9.0 *log( ((mu2>2*LQCD2)?(mu2):(2.0*LQCD2))/LQCD2));
		}
		const double* sigpar;
		double x=0;
	public:
		void set_kinem(const double x){
			this->x=x;
		}
		void init(const double (&par)[]){
			sigpar=par;
		}
		explicit Sigma(void){
		}
		~Sigma(){
		}
		double operator()(const double r,const double x){//,const double Q2,const double*sigpar)const {
		 	this->x=x;
		 	return((*this)(r));
		 }
#if MODEL==0
		double operator()(const double r)const{
			double sigma_0=sigpar[0];
			double lambda=sigpar[1];
			double x_0=sigpar[2];
			if(x_0<1.0e-5||x_0>1.0e-3){
				return 0;
			}
			if(lambda<0.05||lambda>0.95){
				return 0;
			}
			double qs2=pow(x_0/x,lambda); 
#if LAPLACIAN==0
			double val=sigma_0*(1-exp( - pow(r , 2)*qs2/4));
#elif LAPLACIAN==1
			double val=sigma_0*qs2*(1-r*r*qs2/4)*exp(-r*r*qs2/4);
#endif
			return val;
			
			//return( sigma_0*(1-exp( - pow(r * Q0, 2) * pow(x_0/x, lambda)/4)) );	
		}
#elif MODEL==1
	
		//inline PREC alpha(const PREC mu2)const{
		//	return 4.0/(9.0 *log( ((mu2>2*LQCD2)?(mu2):(2.0*LQCD2))/LQCD2));
		//}
		double operator()(const double r){
			double sigma_0=sigpar[0];
			
			//double A_g=sigpar[1];
			//double lambda_g=sigpar[2];
			if(sigpar[1]<0){
				return 0;
			}
			xg.set_xg_parameter(sigpar[1],sigpar[2]);
			double C=sigpar[3];
			double mu02=sigpar[4];
		
			double mu2=C/(r*r)+mu02;
			if(mu2<1){
				return(0);
			}
			double qs2=4*PI*PI*alpha(mu2)*xg(x,mu2)/(3*sigma_0); 
#if LAPLACIAN==0
			double val=sigma_0*(1-exp(-pow(r,2)*qs2/4) );
#elif LAPLACIAN==1
			double val=sigma_0*qs2*(1-r*r*qs2/4)*exp(-r*r*qs2/4);
#endif
			
			if(!std::isfinite(val)){
				printf("%.3e = sigma(%.3e, %.3e;%.3e, %.3e,%.3e, %.3e, %.3e)\n",val,r,x,sigma_0,sigpar[1],sigpar[2],C,mu02);
				printf("%.3e\n",xg(x,mu2));
				getchar();
			} 
			return val;	
		}
#endif		
};

/*class Sigma{
	PREC x=0,Q2=0;
	PREC sigma_0,lambda,x_0;
	std::string key;

	public:	
		int set_kinem(const PREC a,const PREC b){
			if(std::isnan(a+b)+std::isinf(a+b)!=0){
				printf("Sigma:: Q2=%.3le x=%.3le\n",(double)a,(double)b);
				getchar();
			}
			x=a;
			Q2=b;
			return 0;
		}

		explicit Sigma(std::string  type , const PREC(& par)[]){//maybe use struct pointer for parameters.
			if(type=="gbw"){
					sigma_0 =(PREC)par[0];
					lambda	=(PREC)par[1];
					x_0	=(PREC)par[2];
					key=type;
					//printf("sigma_0 %.3le  lambda %.3le x_0 %.3le\n",(double)sigma_0,(double)lambda,(double)x_0);
				}else{
					std::cout<<"unknown model: "<<type<<std::endl;
			}
		}
		~Sigma(){
		}

		PREC sigma(const PREC r)const{
			PREC val=0;
			if(key=="gbw"){
				val=sigma_gbw(r,x,Q2);
			}else{
				
			}
			if(std::isnan(val)+std::isinf(val)!=0){
				printf("sigma_gbw: %.3le encountered \n",(double)val);
				printf("sigma_0 %.3le  lambda %.3le x_0 %.3le\n",(double)sigma_0,(double)lambda,(double)x_0);
				printf("Q2 %.3le  x %.3le \n",(double)Q2,(double)x);
				getchar();
				return 0;
			}else{
				return val;
			}
		}

	private:

		PREC sigma_gbw(const PREC r,const PREC x,const PREC q2)const{
			if(x_0<0){//to avoid nan since migrad might give negative x0...
				return 0;
			}
			PREC result= sigma_0*(1-exp( - pow(r * Q0, 2) * pow(x_0/x, lambda)/4)) ;
			return(result);	

		}
	
};*/
/////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////
/*#if MODEL==0
	typedef Sigma SIGMA;
#elif MODEL==1
	typedef  Laplacian_Sigma SIGMA;
#endif
class Integrand_r{
	SIGMA *sigma_ptr;

	public:
		explicit Integrand_r(const PREC x,const PREC Q2,const PREC mf2,SIGMA&  sig){
			sigma_ptr=&sig;
			set_kinem(x,Q2,mf2);
			//printf("int created %.3e\n",(double)integrand_r(0.1,0.1));
		}
		int set_kinem(const PREC a,const PREC b,const PREC c){
			if(std::isnan(a+b+c)+std::isinf(a+b+c)!=0){
				printf("Integrand:: Q2=%.3le x=%.3le mf2=%.3le\n",(double)a,(double)b,(double)c);
				getchar();
			}
			x=modx(a,b,c);
			Q2=b;
			mf2=c;
#if MODEL==0
			sigma_ptr->set_kinem(x);
#elif MODEL==1
			sigma_ptr->approximate(x);
#endif

			return 0;
		}
	private:
		PREC x, Q2, mf2;
	
	
//	PREC sigma (const PREC r) const {
//		return(sigma_ptr->sigma(r));
//	}
	PREC psisq_f (const PREC z,const PREC r)const  {
		PREC	value;
		PREC	z_bar =  z*z+(1-z)*(1-z);
		PREC	Qsq_bar =  z*(1-z)*Q2+mf2;
		PREC	Qsq2 =  sqrt(Qsq_bar)*r;
		//pow(r,2) is to suppress singularity at r=0, it is compensated by the sigma
		if(Qsq2<1.0e-5){//small er approximation
			value =   (z_bar + ( mf2+ pow(2*z*(1-z),2)* Q2 )*pow(r* log(Qsq2),2) );
		}else{
			PREC	bessel_k0_2 = pow(std::cyl_bessel_k(0,Qsq2),2);
			PREC	bessel_k1_2 = pow(std::cyl_bessel_k(1,Qsq2),2);
			value = pow(r,2) * (z_bar * Qsq_bar * bessel_k1_2 + ( mf2 + pow(2*z*(1-z),2)* Q2 ) * bessel_k0_2);
		}
		PREC result=(3*value)/(2*PI*PI);
		if(std::isnan(result)+std::isinf(result)!=0){
 		 	printf("psiisq_f: %.3le encountered \n",(double)result);
 		 	printf("z %.3le  r %.3le Q2 %.3le mass2 %.3le\n",(double)z,(double)r,(double)Q2,(double)mf2);
 		 	printf("Qsq2 = %.3le, Qsq_bar = %.3le, z_bar = %.3le, value =%.3le   ",(double)Qsq2, (double)Qsq_bar,(double) z_bar,(double)value);
 		 	getchar();
	 	 	return 0;
 		}
		return(result);	
	}
	public:
	//PREC integrand_r(PREC z,PREC r)const{
	PREC operator()(PREC z,PREC r)const{
		PREC jacr=0;
		change_var(r,jacr,R_MIN,R_MAX, 1+Q2);
		PREC jacz=0;
		change_var(z,jacz,0,0.5,10);
		PREC val=(*sigma_ptr)(r)* psisq_f (z, r)/r;
		return(jacr*jacz*2*val);//r^2 comes from photon wave function. just extracted... 2 pi r is angular integration 
	}
};

int F2_integrand_B(const int *__restrict ndim, const PREC *__restrict intv,const int *__restrict ncomp,PREC*__restrict  f, void* __restrict p){
	Integrand_r *integrand=(Integrand_r*)p;
	PREC z=intv[0];
	PREC r=intv[1];
	
	
	PREC res=0;
	
	res+=(2.0/3.0)*integrand[0](z,r);
	res+=(4.0/9.0)*integrand[1](z,r);
	res+=(1.0/9.0)*integrand[2](z,r);

	*f=res;
	return(0);
}*/
