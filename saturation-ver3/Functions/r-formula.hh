////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//   GBW
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Sigma{
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
	
};
/////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////
class Integrand_r{
	Sigma *sigma_ptr;
	
	public:
		explicit Integrand_r(const std::string model,const PREC x,const PREC Q2,const PREC mf2,Sigma&  sig){
			sigma_ptr=&sig;
			set_kinem(x,Q2,mf2);
			//printf("int created %.3e\n",(double)integrand_r(0.1,0.1));
		};
	private:
		PREC x, Q2, mf2;
		int set_kinem(const PREC a,const PREC b,const PREC c){
			if(std::isnan(a+b+c)+std::isinf(a+b+c)!=0){
				printf("Integrand:: Q2=%.3le x=%.3le mf2=%.3le\n",(double)a,(double)b,(double)c);
				getchar();
			}
			x=a;
			Q2=b;
			mf2=c;
			sigma_ptr->set_kinem(x,Q2);
			return 0;
		}
		
	
	PREC sigma (const PREC r) const {
		return(sigma_ptr->sigma(r));
	}
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
		PREC val=sigma(r)* psisq_f (z, r)/r;
		
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
}
