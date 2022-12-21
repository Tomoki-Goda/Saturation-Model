#include<iostream>
#include<cuba.h>
#include<cmath>
#include"./control.h"
#include"./control-default.h"
#include"./constants.h"
extern double INT_PREC;
//#include"./Photon.hh"



#ifndef ALPHA_RUN
	#define ALPHA_RUN 0 
#endif
#ifndef MODX
	#define MODX 0
#endif
#ifndef REUSE
	#define REUSE 0 
#endif

inline double modx(const double x, const double Q2,const  double mf2){
#if MODX==1
	return( (x*(1+4*mf2/Q2)));
#else 
	return( x);
#endif
}

int change_var(double* var,double* jac,const double min, const double max){
	double scale=(max-min);
	double a1=*var;
	double a2=1-a1;
	double den=scale*a2+a1;
	*jac=pow(scale/den,2);
	*var=(min*(scale)*a2+max*a1)/den;
	return 0;
}

class Gluon{
	double sigma_0,lambda,x_0;
	double Q2;
	std::string key;
	
	public:
		Gluon(std::string type,const double *par){
			if(type=="gbw"){
					sigma_0 =par[0];
					lambda	=par[1];
					x_0	=par[2];
					key=type;
			}else{
					std::cout<<"unknown model: "<<type<<std::endl;
			}
		
		}
		~Gluon(){
		}
		int set_kinem(const double b){
			//x=a;
			Q2=b;
			return 0;
		}
	private:
		inline double alpha(const double mu2){
			//const double b0=((11.0*3.0-3.0*2.0)/12.0);
			//const double val=1.0/(b0 *log( ((mu2>2*LQCD2)?(mu2):(2*LQCD2))/LQCD2));
			return 4.0/(9.0 *log( ((mu2>2*LQCD2)?(mu2):(2*LQCD2))/LQCD2));
		}

	public:
		double aF(const double x,const double k2){
			//double lambda=par[1];
			//double x0=par[2];

			double Qs2=pow(x_0/x,lambda);
			double val=3.0/(4*PI*PI)*k2/Qs2*exp(-k2/Qs2);



			if(std::isnan(val)==1){
				return(0);
			}
			return (sigma_0*val) ;
		}
};

class Integrand_kt{
		Gluon * gluptr;
		double x,Q2,mf2;
		
	public:
		Integrand_kt(const std::string type,const double a,const double b,const double c,Gluon* gluon){
			gluptr=gluon;
			set_kinem(a,b,c);
			gluon->set_kinem(b);
		}
		~Integrand_kt(){
		}
		
		int set_kinem(const double a,const double b,const double c){
			x=a;
			Q2=b;
			mf2=c;
			return 0;
		}
		
	private:
		int I_array(const double beta,const  double kappa_t_prime2,const double kt2, double *I)const{
			const double N1=beta*(1-beta)*Q2+mf2;
			const double N2=kappa_t_prime2+pow(1-beta,2)*kt2;
			const double N3=kappa_t_prime2-pow(1-beta,2)*kt2;
			const double N4=kappa_t_prime2+beta*(1-beta)*kt2;
			
			const double nsqrt=sqrt(N1*N1+2*N1*N2+N3*N3);
			I[0]=(N1*N2+N3*N3)/pow(nsqrt,3.0);
			I[1]=(N3-(1-2*beta)*N1)/((N1+N4)*nsqrt);	
			I[2]=(N1+N2)/pow(nsqrt,3.0);
			I[3]=(2*(1-beta))/((N1+N4)*nsqrt);	
			return 0;
		}



		//static 
		inline double inv_z(const double beta, const double kappa_t_prime2,const double kt2)const{
			return( 1+(kappa_t_prime2+mf2)/(beta*(1-beta)*Q2)+kt2/Q2);
		}

		

	public:
		double  integrand(double kt2, double kappa_t_prime2, double beta)const{
			if((1-x)/x *Q2-4*mf2<=0){
				return(0);
			}
			double jac1,jac2,jac3;
			
			const double ktmax=(1-x)/x *Q2-4*mf2;
			change_var(&kt2, &jac1,0,ktmax);
			//double kprimemax=(1-x)/x *Q2/4-mf2;
			const double kprimemax=ktmax/4 ;//just because it is...
			change_var(&kappa_t_prime2, &jac2,0,kprimemax);
			double betamin=sqrt(1-4*(mf2/((1-x)/x*Q2) )), betamax ;
			betamax=(1+betamin)/2;
			betamin=(1-betamin)/2;
			change_var(&beta, &jac3,betamin,betamax);
			if(std::isnan(jac1*jac2*jac3)==1){
				printf("Q2= %.3e, x=%.3e, mf2=%.3e sqrt : %.3e jac1 %.3e jac2 %.3e jac3 %.3e  \n",Q2,x,mf2,1-4*(mf2/((1-x)/x*Q2) ) ,jac1, jac2,jac3);
			}
			const double xz=x*inv_z(beta,kappa_t_prime2,kt2) ;
			if(1-xz<0){
				return 0;		
			}
			
			double val=0;
			double I[4];
			I_array(beta, kappa_t_prime2, kt2, I);

			val+=(beta*beta+pow(1-beta,2))*(I[0]-I[1]);
			val+=(mf2+4*Q2*beta*beta*pow(1-beta,2) )*(I[2]-I[3]);
			val*=gluptr->aF(xz, kt2 );
			
#if ALPHA_RUN==1
			val*=alpha(kappa_t_prime2+kt2+mf2+1)/0.2;
			//printf("%.5e\n", val);
#endif
			return(Q2/(2*PI) *jac1*jac2*jac3* val/kt2 );	
			//return( jac1*jac2*jac3* F2_integrand(beta,kappa_t_prime2,kt2,Q2,mf2,x,par )/kt2 );
		}
};


int F2_integrand_A(const int *ndim, const double* intv,const int *ncomp,double* f, void* p){
	Integrand_kt* param=(Integrand_kt*)p;
	double beta=intv[0];
	double kappa_t_prime2=intv[1];
	double kt2=intv[2];
	double val=0;

	val+= (2.0/3.0)*param[0].integrand( kt2, kappa_t_prime2, beta);
	val+= (4.0/9.0)*param[1].integrand( kt2, kappa_t_prime2, beta);
	val+= (1.0/9.0)*param[2].integrand( kt2, kappa_t_prime2, beta);
	//std::cout<<*f<<std::endl;
	if(std::isnan(val)+std::isinf(val)!=0){
		printf("evaluation failure %.5e\n", val);
	}
	*f=val;
	return 0;
}
/////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////
class Sigma{
	double x=0,Q2=0;
	double sigma_0,lambda,x_0;
	std::string key;

	public:	
		int set_kinem(const double a,const double b){
			x=a;
			Q2=b;
			return 0;
		}

		Sigma(std::string  type , const double* par){//maybe use struct pointer for parameters.
			if(type=="gbw"){
					sigma_0 =par[0];
					lambda	=par[1];
					x_0	=par[2];
					key=type;
				}else{
					std::cout<<"unknown model: "<<type<<std::endl;
			}
		}
		~Sigma(){
		}

		double sigma(const double r)const{
			double val=0;
			if(key=="gbw"){
				val=sigma_gbw(r,x,Q2);
			}else{
				
			}
			if(std::isnan(val)+std::isinf(val)!=0){
				printf("sigma_gbw: %.3e encountered \n",val);
				printf("sigma_0 %.3e  lambda %.3e x_0 %.3e\n",sigma_0,lambda,x_0);
				getchar();
				return 0;
			}else{
				return val;
			}
		}

	private:

		double sigma_gbw(const double r,const double x,const double q2)const{
			if(x_0<0){//to avoid nan since migrad might give negative x0...
				return 0;
			}
			double result= sigma_0*(1-exp( - pow(r * Q0, 2) * pow(x_0/x, lambda)/4)) ;
			return(result);	

		}
	
};
class Integrand_r{
	Sigma *sigma_ptr;
	
	public:
		Integrand_r(std::string model,double x,double Q2,double mf2,Sigma* sig){
			sigma_ptr=sig;
			set_kinem(x,Q2,mf2);
		};
	private:
		double x, Q2, mf2;
		int set_kinem(double a,double b,double c){
			x=a;
			Q2=b;
			mf2=c;
			
			sigma_ptr->set_kinem(x,Q2);
			return 0;
		}
		
	
	double sigma (const double r) const {
		return(sigma_ptr->sigma(r));
	}
	double psisq_f (const double z,const double r)const  {
		double	value;
		double	z_bar =  z*z+(1-z)*(1-z);
		double	Qsq_bar =  z*(1-z)*Q2+mf2;
		double	Qsq2 =  sqrt(Qsq_bar)*r;
		//pow(r,2) is to suppress singularity at r=0, it is compensated by the sigma
		if(Qsq2<1.0e-5){//small er approximation
			value =   (z_bar + ( mf2+ pow(2*z*(1-z),2)* Q2 )*pow(r* log(Qsq2),2) );
		}else{
			double	bessel_k0_2 = pow(std::cyl_bessel_k(0,Qsq2),2);
			double	bessel_k1_2 = pow(std::cyl_bessel_k(1,Qsq2),2);
			value = pow(r,2) * (z_bar * Qsq_bar * bessel_k1_2 + ( mf2 + pow(2*z*(1-z),2)* Q2 ) * bessel_k0_2);
		}
		double result=(3*value)/(2*PI*PI);
		if(std::isnan(result)+std::isinf(result)!=0){
 		 	printf("psiisq_f: %.3e encountered \n",result);
 		 	printf("z %.3e  r %.3e Q2 %.3e mass2 %.3e\n",z,r,Q2,mf2);
 		 	printf("Qsq2 = %.3e, Qsq_bar = %.3e, z_bar = %.3e, value =%.3e   ",Qsq2, Qsq_bar, z_bar,value);
 		 	getchar();
	 	 	return 0;
 		}
		return(result);	
	}
	public:
	double integrand_r(double z,double r)const{
		return((Q2)/(2*PI*r) *sigma(r)* psisq_f (z, r));//r^2 comes from photon wave function. just extracted... 2 pi r is angular integration 
	}
};

int F2_integrand_B(const int *ndim, const double* intv,const int *ncomp,double* f, void* p){
	//double **param=(double**)p;
	
	Integrand_r *param=(Integrand_r*)p;
	
	double z=intv[0];
	double r=intv[1];
	
	double jac=0;
	change_var(&r,&jac,R_MIN,R_MAX);
	double res=0;
	res+=(2.0/3.0)*param[0].integrand_r(z,r);
	res+=(4.0/9.0)*param[1].integrand_r(z,r);
	res+=(1.0/9.0)*param[2].integrand_r(z,r);

	*f=jac*res;
	return(0);
}


double F2_kt(const double x,const  double Q2,const  double mf2,const double* par){
#if R_FORMULA==1
	const int ndim=2;
	Sigma sigma[]={Sigma("gbw",par) ,Sigma("gbw",par) ,Sigma("gbw",par) };
	Integrand_r integrands[]={Integrand_r("gbw",modx(x,Q2,MASS_L2),Q2,MASS_L2,sigma) ,Integrand_r("gbw",modx(x,Q2,MASS_C2),Q2,MASS_C2,sigma+1) ,Integrand_r("gbw",modx(x,Q2,MASS_B2),Q2,MASS_B2,sigma+2) };
	const int key =13;
#else
	const int ndim=3;
	Gluon gluon[]={Gluon("gbw",par) ,Gluon("gbw",par) ,Gluon("gbw",par) };
	Integrand_kt integrands[]={Integrand_kt("gbw", modx(x,Q2,MASS_L2), Q2,MASS_L2, gluon), Integrand_kt("gbw", modx(x,Q2,MASS_C2), Q2, MASS_C2, gluon+1), Integrand_kt("gbw", modx(x,Q2,MASS_B2), Q2, MASS_B2, gluon+2) };
	const int key =11;
#endif
	const long int mineval=pow(10,ndim), maxeval=1/pow(INT_PREC/10,2);//use llChure if larger than ~1.0e+9
	const long int nstart=1.0e+2,nincrease=1.0e+2;

	const int flag= 0+4*0+8*1+16*0+32*0;
	long long int neval;
	int nregions,fail;
	cubareal integral[3],error[3],prob[3];
	
	char statefile[100]="";
 	double result;
 	
 	llCuhre(ndim, 1,
#if R_FORMULA==1
		&F2_integrand_B,
#else
		&F2_integrand_A,
#endif 
		(void*)integrands, 1,INT_PREC,INT_PREC/10, flag, mineval,maxeval, key,statefile,NULL, &nregions, &neval,  &fail, integral, error, prob
	);

/*	
	llVegas(ndim, 1,
#if R_FORMULA==1
			&F2_integrand_B,
#else
			&F2_integrand_A,
#endif 
			(void*)integrands, 1,INT_PREC,INT_PREC,flag, 0,mineval, maxeval, mineval/10, mineval/10, 1000,50, NULL,NULL, &neval,  &fail, integral, error, prob

	);

	if(fail!=0){
		printf("%d",fail );
	}
*/	
/*	llSuave(ndim,1,
#if R_FORMULA==1	
				&F2_integrand_B,
#else
				&F2_integrand_A,
#endif 
				(void*)integrands, 1, INT_PREC,INT_PREC,flag,0, mineval, maxeval,1000,5,25, NULL, NULL,&nregions, &neval, &fail, integral, error, prob
			);
*/	

//	 if(prob[0]>0.1){
//	 	//flag=1;
// 	 	std::cout<<"1. "<<integral[0]<<", "<< error[0]<<",  "<<prob[0]<<std::endl;
// 	 	std::cout<<"1. "<<nregions<<", "<< neval<<",  "<<fail<<std::endl;
//		 Vegas(ndim, 1,
//#if R_FORMULA==1
//			&F2_integrand_B,
//#else
//			&F2_integrand_A,
//#endif 
//			(void*)param, 1,INT_PREC,INT_PREC/100,flag, 0,	mineval, maxeval, nstart, nincrease, mineval/10, mineval/20, NULL,NULL, &neval,  &fail, integral+1, error+1, prob+1
//
//		);
 //		if(prob[1]>0.1){
 //			//std::cout<<"2. "<<integral[1]<<", "<< error[1]<<",  "<<prob[1]<<std::endl;
 //			Suave(ndim,1,
//#if R_FORMULA==1	
//				&F2_integrand_B,
//#else
//				&F2_integrand_A,
//#endif 
//				(void*)param, 1, INT_PREC,INT_PREC/100,flag,0, mineval, maxeval,1000,5,25, NULL, NULL,&nregions, &neval, &fail, integral+2, error+2, prob+2
//			);
//			
// 			if(prob[2]>0.1){
 //				std::cout<<"Cuhre. "<<integral[0]<<", "<< error[0]<<",  "<<prob[0]<<std::endl;
 //				std::cout<<"Vegas. "<<integral[1]<<", "<< error[1]<<",  "<<prob[1]<<std::endl;
 //				std::cout<<"Suave. "<<integral[2]<<", "<< error[2]<<",  "<<prob[2]<<std::endl;
// 			}else{
 //				result=integral[2];
// 			}
// 		}else{
 //			result=integral[1];
// 		}
// 	}else{
 		result=integral[0];
 //	}
 	
 	 //return 0;
 	 if(std::isnan(result)+std::isinf(result)!=0){
 	 	printf("%.3e encountered \n",result);
 	 	return 0;
 	 }
 	 return result;
}



