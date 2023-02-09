#include<iostream>
#include<fstream>
#define TEST 0

#include<cuba.h>
#include<cmath>
#include <vector>
#include"./control.h"
#include"./control-default.h"
#include"./constants.h"
#include"./r-formula.hh"

#if GLUON_APPROX==1
	#include"./interpolation.hh"
#if HANKEL==1
	typedef Hankel_aF Gluon;
#else
	typedef Approx_aF Gluon ;
#endif
	typedef Laplacian_Sigma SIGMA;
#else
	#include"./gluon-gbw.hh"
	typedef Gluon_GBW Gluon ;
	typedef Sigma SIGMA;

#endif
//#include"./clenshaw.h"


//#include"./dgauss.h"
//#include"./clenshaw-curtis.hh"
extern double  INT_PREC ;
extern int N_APPROX;
//#include"./r-formula.hh"
//int CUBACORES=4;
//#include"./Photon.hh"






double change_var(double & var,double &  jac,const double min, const double max,const double c){//This version is (in theory) regular at max->Inf
	double den=( (c==1)?(1):(c+var*(1-c)) );
	jac= ( (min==0.0)?(c*pow(den,-2)*max):(c*pow(den,-2)*(max-min)) ) ;
	var= ( (min==0.0)?(max*var):((max*var+c*min*(1-var)) ))/den;
	//var= (max*var+min*c*(1-var))/den;
	
#if TEST==1	
	if(var>max) {
		if(fabs((var-max)/max)>1.0e-15){
			printf("value below limit %.3e -> %.3e [%.3e, %.3e] diff %.3e, c=%.3e\n",(1-den)/(1-c),var,min,max,var-max, c);
		}
		var=max;
	}else if(var<min){
		if(fabs((min-var)/min)>1.0e-15){
		printf("value below limit %.3e -> %.3e [%.3e, %.3e] diff %.3e, c=%.3e\n",(1-den)/(1-c),var,min,max,min-var, c);
		}
		var=min;
	}
#endif
	return var;
}



//double change_var(double & var,double &  jac,const double min, const double max){
//	return(change_var(var,jac,min,  max,1));
//}



////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//         Angular integral
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename TYPE>class Integrand_kt_phi{
		std::fstream file;
	public:
		Integrand_kt_phi(TYPE& af ){
			//set_kinem(a,b,c);
			gluptr=&af;
			//gluptr->set_kinem(Q2);
			
#if SCATTER==1
			file.open("home/tomoki/Saturation-Model/saturation-ver3/scatter.txt",std::fstream::app);
#endif
		}
		int set_kinem(const double a, const double b, const double c){
			x=a;
			Q2=b;
			mf2=c;
			betamin=1-4*mf2*x/((1-x)*Q2);
			if(betamin<=0.0){
				betamin=0.5;
				betamax=0.5;
			}else{
				betamax=(1+sqrt(betamin))/2;
				betamin=(1-sqrt(betamin))/2;
			}
			kappamax=((1-x)*Q2)/(x*4)-mf2;
			return 0;
		}
		~Integrand_kt_phi(){
#if SCATTER==1
			file.close();	
#endif	
		}

		double operator()(const double x1, const double x2, const double x3, const double x4){
			double k2=x1,kappa2=x2,beta=x3,phi=x4;
			double jac1,jac2,jac3,jac4;
			double k2max;
			if(betamin>=betamax){
				return 0;
			}
			change_var(beta,jac3,betamin,betamax,1);
			change_var(kappa2,jac2,0,kappamax,1+kappamax/pow(Q2,0.25));

			change_var(phi,jac4,0,PI,1);
			jac4*=2;

			k2max=(1-x)/x*Q2-(kappa2+mf2)/(beta*(1-beta));
			if(k2max<=0.0){
				return 0;
			}
			double val=0;
/*
			change_var(k2,jac1,0,k2max,1+k2max/pow(Q2,0.25));
			val=integrand(kappa2,k2,beta,phi)+integrand(kappa2,k2,1-beta,phi+PI);
			val*=jac1*jac2*jac3*jac4;
*/			
			if(kappa2<k2max){
				k2=1-x1*x1;
				change_var(k2,jac1,0,kappa2,1+kappa2/pow(1+Q2,0.5));
				val+=jac1*integrand(kappa2,k2,beta,phi)+integrand(kappa2,k2,1-beta,phi+PI);
				
				k2=x1*x1;
				//change_var(k2,jac1,kappa2,k2max,1+k2max/pow(kappa2*Q2,0.25));
				change_var(k2,jac1,kappa2,k2max,1+(k2max-kappa2));
				val+=jac1 *integrand(kappa2,k2,beta,phi)+integrand(kappa2,k2,1-beta,phi+PI);
				
				val*=2*x1*jac2*jac3*jac4;
			}else{
				change_var(k2,jac1,0,k2max,1+k2max/pow(1+Q2,0.5));
				val=integrand(kappa2,k2,beta,phi)+integrand(kappa2,k2,1-beta,phi+PI);
				val*=jac1*jac2*jac3*jac4;
			}
			
			

#if TEST==1
			if(k2>k2max||std::isnan(val)||std::isinf(val)){
				printf("x1= %.3e x2= %.3e x3= %.3e x4= %.3e\n",x1,x2,x3,x4 );
				printf("kappa2=%.3e k2=%.3e beta=%.3e phi=%.3e\n",kappa2,k2,beta,phi );
				printf(" %.3e <beta<  %.3e\n",betamin,betamax );
				printf("%.3e <k2< %.3e\n",0.0,k2max);
				//printf("%.3e <cos[phi]\n",cosphimin );
				printf("%.3e <kappa2\n",kappamax );
				getchar();
			}
#endif
			return(val);
		}

	private:
		double mf2=0,Q2=0,x=0;
		double betamin=0,betamax=0,kappamax=0;
		TYPE* gluptr=NULL;

		int Ang(const double kappa2,const double k2,const double beta,const double phi, double &A1,double &A2)const{
			const double A=beta*(1-beta)*Q2+mf2;
			const double aA=kappa2+pow(1-beta,2)*k2+A;
			const double b=2*(1-beta)*sqrt(kappa2*k2);
			const double cA=kappa2+pow(beta,2)*k2+A;
			const double d=2*(beta)*sqrt(kappa2*k2);
			const double e=kappa2-beta*(1-beta)*k2;
			const double f=(1-2*beta)*sqrt(kappa2*k2);
			
			const double den=aA+b*cos(phi);
			A1=pow(den,-2)-2*b/((aA*d+cA*b)*(den));
			//A1*=2;
			A2=-A*pow(den,-2)+pow(den,-1)*(1-2*(b*e-aA*f)/(aA*d+cA*b));
			//A2*=2;
			///printf("%.3e %.3e\n",A1,A2);
			return 0;
		}

		
#if SCATTER==1
		double integrand(double kappa2,double k2,double beta,double phi)  {
			file<<kappa2<<"\t"<<k2<<"\t"<<beta<<"\t"<<phi<<"\t"<<Q2<<"\t"<<x<<"\t"<<mf2<<std::endl;
#else
		double integrand(double kappa2,double k2,double beta,double phi) const{
#endif
			double A1,A2;
			Ang(kappa2,k2,beta,phi,A1,A2);
			double xz=0;
			xz=x*(1+(kappa2+mf2)/(beta*(1-beta)*Q2)+k2/Q2);
			//xz=x*(1+(kappa2+mf2)/((1-beta)*Q2)+(kappa2+k2-2*sqrt(kappa2*k2)*cos(phi)+mf2)/(beta*Q2));
			if(xz>1){
				if(fabs(1-xz)>1.0e-10){
					printf("Theta: x/z=%.3e kappa2=%.3e k2=%.3e beta=%.3e phi=%.3e Q2=%.3e x=%.3e mf2=%.3e\n",xz,kappa2,k2,beta,phi,Q2,x,mf2);
					printf("%.3e\n",xz-1);
					getchar();
				}
				return 0;
			}
			double val=0;
			val+=(pow(beta,2)+pow(1-beta,2))*A2;
			val+=(mf2+4*Q2*pow(beta*(1-beta),2))*A1;
			//val*=2;//for using beta<-> 1-beta symmetry
			double mu2=k2*(1+pow(1-beta,2))+kappa2+mf2+2*(1-beta)*sqrt(kappa2*k2)*cos(phi);
			val*=(*gluptr)(xz,k2,mu2);
			return(val/k2);
		}
	
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  phi integrated 
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename TYPE > class Integrand_kt{
		//Gluon * gluptr=NULL;
		TYPE *gluptr=NULL;
		std::fstream file;
	public:
		double  x=0,Q2=0,mf2=0;
		double  betamin=0,betamax=0, k2max=0,kappamax=0;
		explicit Integrand_kt(TYPE & gluon){
			gluptr=&gluon;
			//set_kinem(a,b,c);
			//gluptr->set_kinem(b);
#if SCATTER==1
			//file.open("home/tomoki/Saturation-Model/saturation-ver3/"+gluon.key+"scatter.txt",std::fstream::out);
			file.open("/home/tomoki/Saturation-Model/saturation-ver3/scatter.txt",std::fstream::app);
#endif
		}
		~Integrand_kt(){
#if SCATTER==1
			file.close();	
#endif	
		}
		
		int set_kinem(const double  a,const double  b,const double  c){
			x=a;
			Q2=b;
			mf2=c;
			//gluptr->set_kinem(b);
			k2max=(1-x)/x *Q2-4*mf2;
			kappamax=k2max/4;//just it is
			betamin=sqrt(1-4*(mf2/((1-x)/x*Q2) )) ;
			betamax=(1+betamin)/2;
			betamin=(1-betamin)/2;
			//gluptr->set_max(k2max);	
			//printf("kt2max=%.3e, x=%.3e, Q2=%.3e\n",k2max,x,Q2);
			return 0;
		}
		
	private:
//		int I_array(const double  beta,const  double  kappa_t_prime2,const double  kt2, double (& I)[])const{
		int I_array(const double  beta,const  double  kappa_t_prime2,const double  kt2, double * I )const{
			const double  N1=beta*(1-beta)*Q2+mf2;
			const double  N2=kappa_t_prime2+pow(1-beta,2)*kt2;
			//const double  N3=kappa_t_prime2-pow(1-beta,2)*kt2;
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
			//change_var(kappa2,jac2,0,kappamax,1+kappamax/pow(Q2,0.25));
			change_var(kappa2,jac2,1.0e-10,kappamax,1+kappamax/pow(1+Q2,0.5));

			k2max=(1-x)/x*Q2-(kappa2+mf2)/(beta*(1-beta));
			if(k2max<=0.0){
				return 0;
			}

			double val=0;
/*
			//change_var(k2,jac1,0,k2max,1+k2max/pow(Q2,0.25));
			change_var(k2,jac1,0,k2max,1+k2max/pow(1+Q2,0.5));
			val=integrand(kappa2,k2,beta);//+integrand(kappa2,k2,1-beta);
			val*=jac1*jac2*jac3;
			return val;
*/
			if(kappa2<k2max){
				k2=1-x1*x1;
				change_var(k2,jac1,1.0e-10,kappa2,1+kappa2/pow(1+Q2,0.25));
				val+=jac1*integrand(kappa2,k2,beta);//+integrand(kappa2,k2,1-beta);
				
				k2=x1*x1;
				//change_var(k2,jac1,kappa2,k2max,1+k2max/pow(kappa2*(1+Q2),0.5));
				//change_var(k2,jac1,kappa2,k2max,1+kappamax/pow(1+Q2,0.5));
				change_var(k2,jac1,kappa2,k2max,1+k2max-kappa2);
				val+=jac1 *integrand(kappa2,k2,beta);//+integrand(kappa2,k2,1-beta);
				
				val*=2*x1*jac2*jac3;
			}else{
				change_var(k2,jac1,1.0e-10,k2max,1+k2max/pow(1+Q2,0.25));
				val=integrand(kappa2,k2,beta);//+integrand(kappa2,k2,1-beta);
				val*=jac1*jac2*jac3;
			}


			return val;

		}
		
#if SCATTER==1
		double   integrand(const double  kappa_t_prime2,const double  kt2,const double  beta){
			
			//std::cout<<kappa_t_prime2<<"\t"<<kt2<<"\t"<<beta<<"\t"<<Q2<<"\t"<<x<<"\t"<<mf2<<std::endl;
#else 
		double   integrand(const double  kappa_t_prime2,const double  kt2,const double  beta)const{
#endif
			
			//printf("integrand2\n");
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
				//getchar();
				return 0;		
			}
			double  val=0;
			double  I[4];
			I_array(beta, kappa_t_prime2, kt2, I);

			val+=(beta*beta+pow(1-beta,2))*(I[0]-I[1]);
			val+=(mf2+4*Q2*beta*beta*pow(1-beta,2) )*(I[2]-I[3]);
			double mu2=kt2+kappa_t_prime2+mf2;
			val*=(*gluptr)(xz,kt2,mu2);
			
		
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
#if SCATTER==1
		file<<std::scientific<<kappa_t_prime2<<"\t"<<kt2<<"\t"<<beta<<"\t"<<Q2<<"\t"<<x<<"\t"<<mf2<<"\t"<<val<<std::endl;
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

#if PHI==1
int F2_integrand_A(const int  *ndim,const  double  *intv,const int  *ncomp,double  * f, void * __restrict p){
		Integrand_kt_phi* integrand=(Integrand_kt_phi*)p;
		const double  beta=intv[0];
		const double  kappa_t_prime2=intv[1];
		const double  kt2=intv[2];
		const double phi=intv[3];
		double  val=0;
		val+= (2.0/3.0)*integrand[0]( kt2, kappa_t_prime2, beta,phi);
		val+= (4.0/9.0)*integrand[1]( kt2, kappa_t_prime2, beta,phi);
		val+= (1.0/9.0)*integrand[2]( kt2, kappa_t_prime2, beta,phi);
		*f=val/(4*PI);
		if(std::isnan(val)+std::isinf(val)!=0){
			printf("integrand %.3e encountered beta:%.3e kappa2 %.3e kt2: %.3e phi: %.3e \n",val,beta,kappa_t_prime2,kt2,phi);
			*f=0;
		}
		return 0;
}
#else
int F2_integrand_A(const int  *ndim,const  double  *intv,const int  *ncomp,double  * f, void * __restrict p){
		Integrand_kt<Gluon>* integrand=(Integrand_kt<Gluon>*)p;
		const double  beta=intv[0];
		const double  kappa_t_prime2=intv[1];
		const double  kt2=intv[2];
		double  val=0;
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
#endif

class Integrand_r{
	SIGMA *sigma_ptr;

	public:
		explicit Integrand_r(SIGMA&  sig){
			sigma_ptr=&sig;
			//set_kinem(x,Q2,mf2);
			//printf("int created %.3e\n",(double)integrand_r(0.1,0.1));
		}
		int set_kinem(const double  x,const double  Q2,const double  mf2){
			if(!std::isfinite(x+Q2+mf2)){
				printf("Integrand:: Q2=%.3le x=%.3le mf2=%.3le\n",Q2,x,mf2);
				getchar();
			}
			this->x=modx(x,Q2,mf2);
			this->Q2=Q2;
			this->mf2=mf2;
			
			sigma_ptr->set_kinem(x);

			return 0;
		}
	private:
		double  x, Q2, mf2;
	
	
//	double  sigma (const double  r) const {
//		return(sigma_ptr->sigma(r));
//	}
	double  psisq_f (const double  z,const double  r)const  {
		double 	value;
		double 	z_bar =  z*z+(1-z)*(1-z);
		double 	Qsq_bar =  z*(1-z)*Q2+mf2;
		double 	Qsq2 =  sqrt(Qsq_bar)*r;
		//pow(r,2) is to suppress singularity at r=0, it is compensated by the sigma
		if(Qsq2<1.0e-5){//small er approximation
			value =   (z_bar + ( mf2+ pow(2*z*(1-z),2)* Q2 )*pow(r* log(Qsq2),2) );
		}else{
			double 	bessel_k0_2 = pow(std::cyl_bessel_k(0,Qsq2),2);
			double 	bessel_k1_2 = pow(std::cyl_bessel_k(1,Qsq2),2);
			value = pow(r,2) * (z_bar * Qsq_bar * bessel_k1_2 + ( mf2 + pow(2*z*(1-z),2)* Q2 ) * bessel_k0_2);
		}
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
	public:
	//double  integrand_r(double  z,double  r)const{
	double  operator()(double  z,  double  r)const{
		double  jacr=0;
		change_var(r,jacr,R_MIN,R_MAX,100);// 1+Q2);
		double  jacz=0;
		change_var(z,jacz,0,0.5,10);
		double  val=(*sigma_ptr)(r)* psisq_f (z, r)/r;
		return(jacr*jacz*2*val);//r^2 comes from photon wave function. just extracted... 2 pi r is angular integration 
	}
};

int F2_integrand_B(const int *__restrict ndim, const double  *__restrict intv,const int *__restrict ncomp,double *__restrict  f, void* __restrict p){
	Integrand_r *integrand=(Integrand_r*)p;
	double  z=intv[0];
	double  r=intv[1];
	
	
	double  res=0;
	
	res+=(2.0/3.0)*integrand[0](z,r);
	res+=(4.0/9.0)*integrand[1](z,r);
	res+=(1.0/9.0)*integrand[2](z,r);

	*f=res;
	return(0);
}
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
class F2_kt{
		//int newpar=1;
	
		
#if R_FORMULA==1
			SIGMA sigma[3]={SIGMA() ,SIGMA() ,SIGMA() };

			Integrand_r integrands[3]={
				Integrand_r(sigma[0]) ,
				Integrand_r(sigma[1]) ,
				Integrand_r(sigma[2])
			};
			const int key =13;
			const int ndim=2;
#else//R_FORMULA
			Gluon gluon;//gluon has no flavour dep.
			
			Integrand_kt<Gluon> integrands[3]={
				Integrand_kt( gluon),
				Integrand_kt( gluon),
				Integrand_kt( gluon)
			};
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
		explicit F2_kt(const  double  *par ){
			//this->par=par;
			//printf(" F2 \n");
#if R_FORMULA==1
#if GLUON_APPROX==0
			sigma[0].init(par);
			sigma[1].init(par);
			sigma[2].init(par);
#elif GLUON_APPROX==1
			
			sigma[0].init(N_APPROX+50,par,'s');
			sigma[1].init(N_APPROX+50,par,'s');
			sigma[2].init(N_APPROX+50,par,'s');
			
#endif
#else//R_FORMULA
      ////////////////////////////////////////////
      ///   Momentum space
      ///////////////////////////////////////////
#if GLUON_APPROX==1
			//if( kt2max<Q2*(1-x)/x){//|| (kt2max/10000)>(Q2*(1-x)/x)  ){//EVALUATE ONLY WHEN RANGE IS TOO DIFFERENT
			gluon.init(N_APPROX+50,N_APPROX+50,N_APPROX+50,par);
			gluon.set_max(kt2max);
			//}
#else
			gluon.init(par);
#endif//GLUON_APPROX==1			
#endif//R_FORMULA
		}

		~F2_kt(){
			//printf(" F2 end \n");
			//getchar();
		}
		
		double operator()(const double  x,const  double  Q2,const  double  mf2){
			static unsigned int count;
	//		std::string type="gbw";
			printf("Start F2\t");
			//getchar();
			integrands[0].set_kinem(x,Q2,MASS_L2);
			printf("L+S");
			integrands[1].set_kinem(x,Q2,MASS_C2);
			printf("+C");
			integrands[2].set_kinem(x,Q2,MASS_B2);
			printf("+B\n");
			//getchar();
			
			printf("%d: Cuhre x=%.2e Q2=%.2e\n",count++, x,Q2);
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
				&F2_integrand_B,
#else
				&F2_integrand_A,
				//&F2_integrand_A1,
#endif  
				(void*)integrands,
				 1,INT_PREC ,INT_PREC /10, flag, mineval,maxeval, key,statefile,NULL, &nregions, &neval,  &fail, integral, error, prob
			);
			printf("\033[1A\033[2K\033[1A\033[2K\r");
			//cubawait(&spin);

			result=Q2/(2*PI) *integral[0];
		 
			 if(std::isnan(result)+std::isinf(result)!=0){
				printf("%.3le encountered \n",(double)result);
				return 0;
			 }
			 //printf("end F2 %.3e\n",result);
			 //getchar();
			 return (result);
		}

};











