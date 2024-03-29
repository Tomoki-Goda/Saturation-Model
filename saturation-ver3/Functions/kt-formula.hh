#include<cmath>
#include<iostream>
#include<fstream>
#include<vector>
#include<chrono>
#include<cuba.h>
#include"control.h"
#include"control-default.h"
#include"constants.h"
#include"gluons.hh"
#include"r-formula.hh"
#include"gluon-integrand.hh"
#include"dipole-gluon.hh"
#include"interpolation-gluon.hh"
//#include"types.hh"
#if GLUON_APPROX==1
typedef Approx_aF aF;
#elif MODEL==0
typedef Gluon_GBW aF;
#endif

///////////////////////////////////////////////////////////////////////////////
//
//
//  Factorization formula of F2 in Kimber's PhD Thesis at Durham
//  phi integrated 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if R_FORMULA==0
template<typename TYPE > class Integrand_kt{
		//Gluon * gluptr=NULL;
		TYPE *gluptr=NULL;
		//std::fstream file;
	public:
		double  x=0,Q2=0,mf2=0;
		double  betamin=0,betamax=0, k2max=0,kappamax=0;
		double 	k2min=KT2_MIN;
		Integrand_kt& operator=(const Integrand_kt& rhs){
			x=rhs.x;
			Q2=rhs.Q2;
			mf2=rhs.mf2;
			betamax=rhs.betamin;
			k2max=rhs.k2max;
			kappamax=rhs.kappamax;
			gluptr=rhs.gluptr;
			return *this;
		}
		
		explicit Integrand_kt(TYPE & gluon){
			gluptr=&gluon;
		}

		~Integrand_kt(){
			//printf("Integrand_kt\n");
		}
		
		int set_kinem(const double  a,const double  b,const double  c){
			x=a;
			Q2=b;
			mf2=c;
			//some restrictions due to the theta function.
			k2max=(1-x)/x *Q2-4*mf2;
			kappamax=k2max/4;//just it is
			betamin=sqrt(1-4*(mf2/((1-x)/x*Q2) )) ;
			betamax=(1+betamin)/2;
			betamin=(1-betamin)/2;
			/*
			if(betamin>betamax){
				printf("beta range error: [%.2e, %.2e] x=%.2e Q2=%.2e mf2=%.2e\n",betamin betamax,x,Q2,mf2);
			}
			if(k2min>k2max){
				printf("k2 range error: [%.2e, %.2e] x=%.2e Q2=%.2e mf2=%.2e\n",k2min k2max,x,Q2,mf2);
			}
			if(k2min>kappamax){
				printf("kappa2 range error: [%.2e, %.2e] x=%.2e Q2=%.2e mf2=%.2e\n",k2min kappamax,x,Q2,mf2);
			}*/
			return 0;
		}
		
	private:
		int I_array(const double  beta,const  double  kappa_t_prime2,const double  kt2, double * I )const{
			const double  N1=beta*(1-beta)*Q2+mf2;
			const double  N2=kappa_t_prime2+pow(1-beta,2)*kt2;
			const double  N3=2*kappa_t_prime2-N2;
			const double  N4=kappa_t_prime2+beta*(1-beta)*kt2;
			const double  nsqrt=sqrt(N1*N1+2*N1*N2+N3*N3);
			const double  den1=pow(nsqrt,3.0);
			const double  den2=(N1+N4)*nsqrt;
			I[0]=(N1*N2+N3*N3)/den1;
			I[1]=(-(1-2*beta)*N1+N3)/den2;	
			I[2]=(N1+N2)/den1;
			I[3]=(2*(1-beta))/den2;	
			return 0;
		}
		inline double  inv_z(const double  beta, const double  kappa_t_prime2,const double  kt2)const{
			return( 1+(kappa_t_prime2+mf2)/(beta*(1-beta)*Q2)+kt2/Q2);
		}
	public:
		inline double  kt2_max(const  double  kappa_t_prime2,const double  beta)const{
			double  val= (1-x)/x *Q2-(kappa_t_prime2+mf2)/(beta*(1-beta)) ;
			return val;
		}
		double operator()(const double x1, const double x2, const double x3){
			double k2=x1,kappa2=x2,beta=x3;
			double jac1,jac2,jac3;
			double k2max;
			if(betamin>=betamax){
				return 0;
			}
			
			change_var(beta,jac3,betamin,betamax,1);
			change_var(kappa2,jac2,k2min,kappamax,1+kappamax/pow(1+Q2,0.5));

			k2max=(1-x)/x*Q2-(kappa2+mf2)/(beta*(1-beta));
			
			if(k2max<=k2min){
				return 0;
			}

			double val=0;
					
			if(kappa2<k2max){
				//cutting space along kappa2=k2.
				// not strictly necesssary but possiblly advantageous.
				// there is a ridge at kappa2=k2 .
				k2=1-x1*x1;
				change_var(k2,jac1,k2min,kappa2,1+kappa2/pow(1+Q2,0.25));
				val+=jac1*integrand(kappa2,k2,beta);//+integrand(kappa2,k2,1-beta);
				
				k2=x1*x1;
				change_var(k2,jac1,kappa2,k2max,1+k2max-kappa2);
				val+=jac1 *integrand(kappa2,k2,beta);//+integrand(kappa2,k2,1-beta);
				
				val*=2*x1*jac2*jac3;
			}else{
				change_var(k2,jac1,k2min,k2max,1+k2max/pow(1+Q2,0.25));
				val=integrand(kappa2,k2,beta);//+integrand(kappa2,k2,1-beta);
				val*=jac1*jac2*jac3;
			}


			return val;

		}
		
		double   integrand(const double  kappa_t_prime2,const double  kt2,const double  beta)const{
			const double  xz=x*inv_z(beta,kappa_t_prime2,kt2) ;
			if(xz>1.0){
#if TEST==1
				if(1-xz<-1.0e-15){
					printf("Thetafunction error \n");
					printf("Q2= %.3le, x=%.3le, mf2=%.3le sqrt : %.3le  \n",
					(double)Q2,(double)x,(double)mf2,(double)(1-4*(mf2/((1-x)/x*Q2) )));
					printf("kt2=%.3e, kappa_t_prime2=%.3e,beta=%.3e  \n",
					(double)kt2,(double)kappa_t_prime2,(double)beta);
					printf("%.3e\n",1-xz);
					getchar();
				}
#endif
				return 0;		
			}
			double  val=0;
			double  I[4];
			I_array(beta, kappa_t_prime2, kt2, I);

			val+=(beta*beta+pow(1-beta,2))*(I[0]-I[1]);
			val+=(mf2+4*Q2*beta*beta*pow(1-beta,2) )*(I[2]-I[3]);
			double mu2=kt2+kappa_t_prime2+mf2;
			val*=(*gluptr)(xz,kt2,mu2);
#if ALPHA_RUN!=0
			val*=alpha(mu2+MU02)/0.2;
#endif
			//val*=(*gluptr)(kt2,mu2);
			
		
			val= val/kt2;
#if TEST==1
			if(std::isnan(val)+std::isinf(val)!=0){
				printf("integrand2 %f encountered\n",val);
				printf("Q2= %.3le, x=%.3le, mf2=%.3le sqrt : %.3le  \n",
				(double)Q2,(double)x,(double)mf2,(double)(1-4*(mf2/((1-x)/x*Q2) )) );
				printf("kt2=%.3e, kappa_t_prime2=%.3e,beta=%.3e  \n",
				(double)kt2,(double)kappa_t_prime2,(double)beta);
				getchar();
			}
#endif
			return(val);
		}
		
		int F2_integrand_A(const int  *ndim,const  double  *intv,const int  *ncomp,double  * f, void * __restrict p){
			Integrand_kt* integrand=(Integrand_kt*)p;
			const double  beta=intv[0];
			const double  kappa_t_prime2=intv[1];
			const double  kt2=intv[2];	
			//double  val=0;
			*f=(*this)( kt2, kappa_t_prime2, beta);
			return 0;
		}

};

template <typename Gluon> int F2_integrand_A(const int  *ndim,const  double  *intv,const int  *ncomp,double  * f, void * __restrict p){
		Integrand_kt<Gluon>* integrand=(Integrand_kt<Gluon>*)p;
		const double  beta=intv[0];
		const double  kappa_t_prime2=intv[1];
		const double  kt2=intv[2];
		double  val=0;
		//sum over flavours
		val+= (2.0/3.0)*integrand[0]( kt2, kappa_t_prime2, beta);
		val+= (4.0/9.0)*integrand[1]( kt2, kappa_t_prime2, beta);
		val+= (1.0/9.0)*integrand[2]( kt2, kappa_t_prime2, beta);
		*f=val;
		if(std::isnan(val)+std::isinf(val)!=0){
			printf("integrand %.3e encountered beta:%.3e kappa2 %.3e kt2: %.3e \n",val,beta,kappa_t_prime2,kt2);
			*f=0;
		}
		return 0;
}

////////////////////////////////////////////////////////////////////
//
//
//  Dipole factorization formula (See Golec-Biernat Wusthoff etc)
//
//
/////////////////////////////////////////////////////////////////////

#elif R_FORMULA==1
template <typename Sig> class Integrand_r{
	Sig *sigma_ptr;
	public:
		explicit Integrand_r(Sig&  sig){
			sigma_ptr=&sig;
		}
		int set_kinem(const double x,const double Q2,const double  mf2){
		//	if(!std::isfinite(x+Q2+mf2)){
		//		printf("Integrand:: Q2=%.3le x=%.3le mf2=%.3le\n",Q2,x,mf2);
		//		getchar();
		//	}
			this->x=modx(x,Q2,mf2);
			//this->x=x;
			if((this->x)>1){
				printf("set_kinem:: x too large %.3e\n",x);
			}
			this->Q2=Q2;
			this->mf2=mf2;
			//if(typeof(*sigma_ptr)==typeof(Sigma_BGK)){
				sigma_ptr->set_x(this->x);
				//sigma_ptr->set_x(x);
			//}
			return 0;
		}
		double  operator()(double  z,  double  r)const{
			double  jacr=0;
			change_var(r,jacr,R_MIN,R_MAX,100);// 1+Q2);
			double  jacz=0;
			change_var(z,jacz,0,0.5,10);
			double  val=(*sigma_ptr)(x,r)* psisq_f(z, r)/r;
			return(jacr*jacz*2*val);//r^2 comes from photon wave function. just extracted... 2 pi r is angular integration 
		}
	private:
		double  x,Q2, mf2;
		double  psisq_f (const double  z,const double  r)const  {
			double 	value;
			double 	z_bar =  z*z+(1-z)*(1-z);
			double 	Qsq_bar =  z*(1-z)*Q2+mf2;
			double 	Qsq2 =  sqrt(Qsq_bar)*r;
			//pow(r,2) is to suppress singularity at r=0, it is compensated by the sigma
//			if(Qsq2<1.0e-5){//small er approximation
//				value =   (z_bar + ( mf2+ pow(2*z*(1-z),2)* Q2 )*pow(r* log(Qsq2),2) );
//			}else{

				double 	bessel_k0_2 = pow(std::cyl_bessel_k(0,Qsq2),2);
				double 	bessel_k1_2 = pow(std::cyl_bessel_k(1,Qsq2),2);
				value = pow(r,2) * (z_bar * Qsq_bar * bessel_k1_2 + ( mf2 + pow(2*z*(1-z),2)* Q2 ) * bessel_k0_2);
//				if(Qsq2<2.5e-5){
//					if(fabs(value - (z_bar + ( mf2+ pow(2*z*(1-z),2)* Q2 )*pow(r* log(Qsq2),2) ))>1.0e-10*value ){
//						printf("psisq_f discontinuous\n%2e %.2e\n",value, (z_bar + ( mf2+ pow(2*z*(1-z),2)* Q2 )*pow(r* log(Qsq2),2)) );
//					}
//				}
//			}

			double  result=(3*value)/(2*PI*PI);
			if(std::isnan(result)+std::isinf(result)!=0){
 			 	printf("psiisq_f: %.3le encountered \n",(double)result);
 			 	printf("z %.3le  r %.3le Q2 %.3le mass2 %.3le\n",(double)z,(double)r,(double)Q2,(double)mf2);
 			 	printf("Qsq2 = %.3le, Qsq_bar = %.3le, z_bar = %.3le, value =%.3le   ",(double)Qsq2, (double)Qsq_bar,(double) z_bar,(double)value);
 			 	getchar();
		 	 	return 0;
 			}
			return(result);	
		}
};

template <typename Sig>int F2_integrand_B(const int *__restrict ndim, const double  *__restrict intv,const int *__restrict ncomp,double *__restrict  f, void* __restrict p){
	Integrand_r<Sig> *integrand=(Integrand_r<Sig>*)p;
	double  z=intv[0];
	double  r=intv[1];
	double  res=0;
	//sum over flavours
	res+=(2.0/3.0)*integrand[0](z,r);
	res+=(4.0/9.0)*integrand[1](z,r);
	res+=(1.0/9.0)*integrand[2](z,r);

	*f=res;
	return(0);
}
#endif


////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void llTest(const int ndim, const int ncomp,
  integrand_t integrand, void *userdata, const long long int nvec,
  const cubareal epsrel, const cubareal epsabs,
  const int flags,
  const long long int mineval, const long long int maxeval,
  const int key,
  const char *statefile, void *spin,
  int *nregions, long long int *neval, int *fail,
  cubareal integral[], cubareal error[], cubareal prob[]){
  printf("calling\n");
  
  }
  
///////////////////////////////////////////////////////
//
// Integration is pertformed by the 'Cuba' library
//
//////////////////////////////////////////////////////
template <typename T> class F2_kt{
			T *integrands;
#if R_FORMULA==1
			const int key =13;
			const int ndim=2;
#else//R_FORMULA
			const int key =11;
			const int ndim=3;
			const double kt2max=5.0e+4;
			//printf(" %.3e \n",Q2*(1-x)/x);	
#endif//R_FORMULA	
      ////////////////////////////////////////////
      ///   Position space
      ///////////////////////////////////////////
			//const double*__restricted par;
	public: 
		//explicit F2_kt(const F2_kt & init){
		//	this->par=init.par;
		//	this->integrands=init.integrands;
		//	for(int i=0;i<3;i++){
		//		(this->integrands)[i]=(init.integrands)[i];
		//	}
		//}
		
		explicit F2_kt(T* integrands ){
			this->integrands=integrands;
		}

		~F2_kt(){
		}
		
		double operator()(const double  x,const  double  Q2,const  double  mf2){
			static unsigned int count;
			//printf("Start F2\t");
			//getchar();
			integrands[0].set_kinem(x,Q2,MASS_L2);
			//printf("L+S");
			integrands[1].set_kinem(x,Q2,MASS_C2);
			//printf("+C");
			integrands[2].set_kinem(x,Q2,MASS_B2);
			//printf("+B\n");
			//getchar();
			
			//printf("%d: Cuhre x=%.2e Q2=%.2e\n",count++, x,Q2);
			const long long int mineval=pow(15,ndim), maxeval=1/pow(INT_PREC /10,2);//use llChure if larger than ~1.0e+9
			const long long int nstart=1.0e+2,nincrease=1.0e+2;
			long long int neval=0;
			//const int mineval=pow(10,ndim), maxeval=1/pow(INT_PREC /10,2);//use llChure if larger than ~1.0e+9
			//const int nstart=1.0e+2,nincrease=1.0e+2;
			//int neval=0;

			const int flag= 0+4*0+8*1+16*0+32*0;
			
			int nregions=0,fail=0;
			//cubareal
			double  integral[3]={0},error[3]={0},prob[3]={0};
			//int spin=0;
			char statefile[100]="";
			double  result=0;
			llCuhre(ndim, 1,
			//llTest(ndim, 1,
#if R_FORMULA==1
//	#if SIGMA_APPROX==-2||SIGMA_APPROX==1
//				&F2_integrand_B<DSIGMA>,
//	#else
				&F2_integrand_B<SIGMA>,
//	#endif
#else

				&F2_integrand_A<aF>,
#endif  
				(void*)integrands,
				 1,INT_PREC ,INT_PREC /10, flag, mineval,maxeval, key,NULL,NULL, &nregions, &neval,  &fail, integral, error, prob
			);
			//printf("\033[1A\033[2K\033[1A\033[2K\r");
			//printf("\033[1A\033[2K\r");
			//cubawait(&spin);

			result=Q2/(2*PI) *integral[0];
		 
			 if(!std::isfinite(result)){
				printf("%.3le encountered \n",(double)result);
				return 0;
			 }
			 //printf("end F2 %.3e\n",result);
			 //getchar();
			 return (result);
		}

};











