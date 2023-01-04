#include<iostream>

#include<cuba.h>
#include<cmath>
#include"./control.h"
#include"./control-default.h"
#include"./constants.h"

//#include"./clenshaw.h"
//#include"./dgauss.h"
#include"./clenshaw-curtis.hh"
extern PREC INT_PREC;

//int CUBACORES=4;
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

inline PREC modx(const PREC x, const PREC Q2, const  PREC mf2){
	//if(std::isnan(x+Q2+mf2+mf2/Q2)+std::isinf(x+Q2+mf2+mf2/Q2)!=0){
	//			printf("modx:: Q2=%.3le x=%.3le mf2=%.3le\n",(double)Q2,(double)x,(double)mf2);
	//			getchar();
	//}
#if MODX==1
	return( (x*(1+4*mf2/Q2)));
#else 
	return( x);
#endif
}

double change_var(PREC& var,PREC&  jac,const PREC min, const PREC max){//This version is regular at max->Inf
	PREC scale=(max-min);
	jac=scale*2*var;
	var=min+scale*pow(var,2);
	return var;
}
double change_var2(PREC & var,PREC &  jac,const PREC min, const PREC max){//This version is regular at max->Inf
	PREC scale=(max-min);
	PREC den=scale*(1-var)+var;
	jac=pow(scale/den,2);
	var=(min*(scale)*(1-var)+max*var)/den;	
	return var;
}
/*PREC change_var(PREC* __restrict var,PREC*  __restrict jac,const PREC min, const PREC max){
			PREC scale=(max-min);
			*jac=scale;
			*var=min+scale*(*var);
			return *var;
}*/

class Gluon{
	PREC sigma_0=0,lambda=0,x_0=0;
	PREC Q2=0;
	std::string key;
	
	public:
		explicit Gluon(std::string &type,const PREC (&par)[]){
				printf("gluon created\n");
				key=type;	
				set_par(par);	
		}
		void set_par(const PREC(&par)[]){
			if(key=="gbw"){
				sigma_0 =(PREC)par[0];
				lambda	=(PREC)par[1];
				x_0	=(PREC)par[2];
			}else{
				std::cout<<"unknown model: "<<key<<std::endl;
			}

		}
		~Gluon(){
		}
		int set_kinem(const PREC b){
			//x=a;
			Q2=b;
			return 0;
		}
	public:
		inline PREC alpha(const PREC mu2){
			//const PREC b0=((11.0*3.0-3.0*2.0)/12.0);
			//const PREC val=1.0/(b0 *log( ((mu2>2*LQCD2)?(mu2):(2*LQCD2))/LQCD2));
			return 4.0/(9.0 *log( ((mu2>2*LQCD2)?(mu2):(2.0*LQCD2))/LQCD2));
		}

	public:
		PREC aF(const PREC x,const PREC k2){
			if(x_0<0||x_0>0.01){
				return 0;
			}
			if(lambda<0||lambda>2){
				return 0;
			}
			PREC Qs2=pow(x_0/x,lambda);
			PREC val=3.0/(4*PI*PI)*k2/Qs2*exp(-k2/Qs2);
			if(std::isnan(val)==1){
				return(0);
			}
			return (sigma_0*val) ;
		}
};

class Integrand_kt{
		Gluon * gluptr=NULL;
		//PREC x,Q2,mf2;
		//PREC jackt2 jackappat2,jacbeta;
		
	public:
		PREC x=0,Q2=0,mf2=0;
		PREC betamin=0,betamax=0, ktmax=0,kprimemax=0;
		explicit Integrand_kt(const std::string &type,const PREC a,const PREC b,const PREC c,Gluon & gluon){
			gluptr=&gluon;
			set_kinem(a,b,c);
		}
		void set(const PREC a,const PREC b,const PREC c){
			//gluptr=&gluon;
			set_kinem(a,b,c);
			//gluon.set_kinem(b);
		}
		~Integrand_kt(){
		}
		
		int set_kinem(const PREC a,const PREC b,const PREC c){
			x=a;
			Q2=b;
			mf2=c;
			
			gluptr->set_kinem(b);
			ktmax=(1-x)/x *Q2-4*mf2;
			//change_var(&kt2, &jackt2,0,ktmax);
			
			kprimemax=(1-x)/x *Q2/4-mf2;
			//change_var(&kappa_t_prime, &jackappa,0,kprimemax);
			
			betamin=sqrt(1-4*(mf2/((1-x)/x*Q2) )) ;
			betamax=(1+betamin)/2;
			betamin=(1-betamin)/2;
			//betamin=1.0/2.0;
			//change_var(&beta, &jacbeta,betamin,betamax);
			
			return 0;
		}
		
	private:
		int I_array(const PREC beta,const  PREC kappa_t_prime2,const PREC kt2, PREC(& I)[])const{
			const PREC N1=beta*(1-beta)*Q2+mf2;
			const PREC N2=kappa_t_prime2+pow(1-beta,2)*kt2;
			const PREC N3=kappa_t_prime2-pow(1-beta,2)*kt2;
			const PREC N4=kappa_t_prime2+beta*(1-beta)*kt2;
			
			const PREC nsqrt=sqrt(N1*N1+2*N1*N2+N3*N3);
			I[0]=(N1*N2+N3*N3)/pow(nsqrt,3.0);
			I[1]=(N3-(1-2*beta)*N1)/((N1+N4)*nsqrt);	
			I[2]=(N1+N2)/pow(nsqrt,3.0);
			I[3]=(2*(1-beta))/((N1+N4)*nsqrt);	
			return 0;
		}
		inline PREC inv_z(const PREC beta, const PREC kappa_t_prime2,const PREC kt2)const{
			return( 1+(kappa_t_prime2+mf2)/(beta*(1-beta)*Q2)+kt2/Q2);
		}
	public:
		inline PREC kt2_max(const  PREC kappa_t_prime2,const PREC beta)const{
			PREC val= (1-x)/x *Q2-(kappa_t_prime2+mf2)/(beta*(1-beta)) ;
			return val;
		}
		PREC  integrand(const PREC x1,const  PREC x2,const PREC x3)const{
			//printf("integrand\n");
			PREC kt2=x1, kappa_t_prime2=x2, beta=x3;
			PREC jac1=1,jac2=1,jac3=1;
			
			change_var2(kappa_t_prime2, jac2,1.0e-15,kprimemax);
			change_var(beta, jac3,betamin,betamax);
			
			PREC ktmax=kt2_max(kappa_t_prime2,  beta);
			if(1.0e-15>=ktmax){
				return 0;
			}
			change_var2(kt2, jac1,1.0e-15,ktmax);
			/*if((std::isnan(jac1*jac2*jac3*kappa_t_prime2*beta*kt2 )
			+std::isinf(jac1*jac2*jac3*kappa_t_prime2*beta*kt2))!=0){
				printf("integrand:: %f encountered",jac1*jac2*jac3*kappa_t_prime2*beta*kt2);
				printf("Q2= %.3le, x=%.3le, mf2=%.3le sqrt : %.3le jac1 %.3le jac2 %.3le jac3 %.3le  \n",
				(double)Q2,(double)x,(double)mf2,(double)(1-4*(mf2/((1-x)/x*Q2) )),
				(double)jac1,(double) jac2,(double)jac3);
				printf("x1= %.3le, x2=%.3le, x3=%.3le \n",(double)jac1,(double) jac2,(double)jac3);
				getchar();
			}*/
			if(kt2<0||kt2>ktmax){
				printf("wrong %.3e -> %.3e <ktmax=%.3e\n",x1, kt2,ktmax );
			}
			
			PREC val;
			val=integrand2( beta,kappa_t_prime2, kt2);
			return(jac1*jac2*jac3*val);
		}
		PREC  integrand01(const PREC x1,const  PREC x2,const PREC x3)const{
			PREC kt2=x1, kappa_t_prime2=x2, beta=x3;
			PREC jac1=1,jac2=1,jac3=1;
			change_var2(kappa_t_prime2, jac2,1.0e-15,kprimemax);
			change_var(beta, jac3,betamin,betamax);
			PREC ktmax=kt2_max(kappa_t_prime2,  beta);
			if(1.0e-15>=kappa_t_prime2||(ktmax<kappa_t_prime2)){
				return 0;
			}
			change_var2(kt2, jac1,1.0e-15,kappa_t_prime2);
			if(kt2<0||kt2>ktmax){
				printf("wrong %.3e -> %.3e <ktmax=%.3e\n",x1, kt2,ktmax );
			}
			
			PREC val;
			val=integrand2( beta,kappa_t_prime2, kt2);
			return(jac1*jac2*jac3*val);
		}
		PREC integrand02(const PREC x1,const  PREC x2,const PREC x3)const{
			PREC kt2=x1, kappa_t_prime2=x2, beta=x3;
			PREC jac1=1,jac2=1,jac3=1;
			
			change_var2(kappa_t_prime2, jac2,1.0e-15,kprimemax);
			change_var(beta, jac3,betamin,betamax);
			
			PREC ktmax=kt2_max(kappa_t_prime2,  beta);
			if(kappa_t_prime2>=ktmax||kappa_t_prime2>ktmax){
				return 0;
			}
			change_var2(kt2, jac1,kappa_t_prime2,ktmax);
			if(kt2<0||kt2>ktmax){
				printf("wrong %.3e -> %.3e <ktmax=%.3e\n",x1, kt2,ktmax );
			}
			
			PREC val;
			val=integrand2( beta,kappa_t_prime2, kt2);
			return(jac1*jac2*jac3*val);
		}
		
		PREC  integrand2(const PREC beta, const PREC kappa_t_prime2,const PREC kt2)const{
			
			//printf("integrand2\n");
			const PREC xz=x*inv_z(beta,kappa_t_prime2,kt2) ;
			if(xz>1.0){
				if(1-xz<-1.0e-15){
					printf("Thetafunction error \n");
					printf("Q2= %.3le, x=%.3le, mf2=%.3le sqrt : %.3le  \n",
					(double)Q2,(double)x,(double)mf2,(double)(1-4*(mf2/((1-x)/x*Q2) )));
					printf("kt2=%.3e, kappa_t_prime2=%.3e,beta=%.3e  \n",
					(double)kt2,(double)kappa_t_prime2,(double)beta);
					printf("%.3e\n",1-xz);
					getchar();
				}
				//getchar();
				return 0;		
			}
			PREC val=0;
			PREC I[4];
			I_array(beta, kappa_t_prime2, kt2, I);

			val+=(beta*beta+pow(1-beta,2))*(I[0]-I[1]);
			val+=(mf2+4*Q2*beta*beta*pow(1-beta,2) )*(I[2]-I[3]);
			val*=gluptr->aF(xz, kt2 );
			
#if ALPHA_RUN==1
			val*=gluptr->alpha(kappa_t_prime2+kt2+mf2+1)/0.2;
			//printf("%.5le\n", val);
#endif			
			val=Q2/(2*PI) * val/kt2;
			if(std::isnan(val)+std::isinf(val)!=0){
				printf("integrand2 %f encountered\n",val);
				printf("Q2= %.3le, x=%.3le, mf2=%.3le sqrt : %.3le  \n",
				(double)Q2,(double)x,(double)mf2,(double)(1-4*(mf2/((1-x)/x*Q2) )) );
				printf("kt2=%.3e, kappa_t_prime2=%.3e,beta=%.3e  \n",
				(double)kt2,(double)kappa_t_prime2,(double)beta);
				getchar();
			}
			return(val);
		}
			
		PREC  integrand_secdec(const PREC x1,const  PREC x2,const PREC x3)const{
			double jacb=2*x3;
			double beta=pow(x3,2);
			
			double val=0, x12,kappa_t_prime2,kt2;
			double jac , x122;
			x12=(1-x1*x1)*x2;
			jac=2*x1*x2;
			//x122=pow(x12,2);
			val+= jac*integrand(x12,x2, beta);
			//x12=x1*x2;
			x12=(1-x2*x2)*x1;
			//x122=pow(x12,2);
			jac=2*x1*x2;
			val+=jac*integrand(x1,x12, beta);
			return(jacb*val);
		}
		
};
struct int_param{

	Integrand_kt* int_ptr=NULL;
	const PREC *par=NULL;
	
};
PREC F2_integrand_A0(PREC * __restrict Kt2,void* __restrict p){
	//printf("start A0\n");
	int_param* param=(int_param*)p;
	PREC beta=param->par[0];
	PREC kappa_t_prime2=param->par[1];
	PREC kt2=*Kt2;
	PREC val=0;
	//val+=(param->int_ptr)->integrand_secdec( kt2, kappa_t_prime2, beta);
	val+=(param->int_ptr)->integrand( kt2, kappa_t_prime2, beta);
	if(std::isnan(val)+std::isinf(val)!=0){
		printf("evaluation failure %.5le\n", (double)val);
	}
	//printf("end A0\n");
	return val;
}


int F2_integrand_A1(const int * __restrict ndim, const PREC* __restrict  intv,const int * __restrict ncomp,PREC* __restrict f, void* __restrict p){
	//printf("start A1\n");
	static int count;
	count++;
	PREC kappa2=intv[1],beta=intv[0], dummy3;
	PREC ktmax;
	Clenshaw_Curtis integrator(16);
	integrator.DIV=25;
	integrator.MAX_RECURSION=4;
	//integrator.REV=5;
	Integrand_kt* param=(Integrand_kt*)p;
	int_param param_str;
	param_str.par=intv;
	PREC val=0,val0;
	
	//ktmax=(param_str.int_ptr)->kt2_max(intv[1],intv[0]);
	
	param_str.int_ptr=param;
	val0=integrator.integrate(&F2_integrand_A0,(void*)&param_str,0.0,1.0,(int)(0.1/INT_PREC),INT_PREC/1,INT_PREC/10);
	
	if(integrator.ERROR==1){
		printf("kap=%.3e, b=%.3e,Q2=%.3e, x=%.3e, mf2=%.3e, ktmax=%.3e\n",
		change_var2(kappa2,dummy3,0,(param_str.int_ptr)->kprimemax),change_var(beta,dummy3,(param_str.int_ptr)->betamin,(param_str.int_ptr)->betamax),
		(param_str.int_ptr)->Q2,(param_str.int_ptr)->x,(param_str.int_ptr)->mf2,(param_str.int_ptr)->kt2_max(intv[1],intv[0]));
		getchar();
	}
	val+=(2.0/3.0)*val0;
	
	param_str.int_ptr=param+1;
	val0=integrator.integrate(&F2_integrand_A0,(void*)&param_str,0.0,1.0,(int)(0.1/INT_PREC),INT_PREC/1,INT_PREC/10);
	if(integrator.ERROR==1){
		printf("kap=%.3e, b=%.3e,Q2=%.3e, x=%.3e, mf2=%.3e, ktmax=%.3e\n",
		change_var2(kappa2,dummy3,0,(param_str.int_ptr)->kprimemax),change_var(beta,dummy3,(param_str.int_ptr)->betamin,(param_str.int_ptr)->betamax),
		(param_str.int_ptr)->Q2,(param_str.int_ptr)->x,(param_str.int_ptr)->mf2,(param_str.int_ptr)->kt2_max(intv[1],intv[0]));
		getchar();
	}
	val+=(4.0/9.0)*val0;
	
	param_str.int_ptr=param+2;
	val0=integrator.integrate(&F2_integrand_A0,(void*)&param_str,0.0,1.0,(int)(0.1/INT_PREC),INT_PREC/1,INT_PREC/10);
	if(integrator.ERROR==1){
		printf("kap=%.3e, b=%.3e,Q2=%.3e, x=%.3e, mf2=%.3e, ktmax=%.3e\n",
		change_var2(kappa2,dummy3,0,(param_str.int_ptr)->kprimemax),change_var(beta,dummy3,(param_str.int_ptr)->betamin,(param_str.int_ptr)->betamax),
		(param_str.int_ptr)->Q2,(param_str.int_ptr)->x,(param_str.int_ptr)->mf2,(param_str.int_ptr)->kt2_max(intv[1],intv[0]));
		//printf("kap=%.3e, b=%.3e,Q2=%.3e, x=%.3e, mf2=%.3e\n  ",intv[1],intv[0],param->Q2,param->x,param->mf2);
		getchar();
	}
	if(std::isnan(val)+std::isinf(val)!=0){
		printf("evaluation failure %.5le\n", (double)val);
	}
	val+=(1.0/9.0)*val0;
	
	//printf("%d val=%.2e\n",count,val);
	*f=val;
	return 0;
}
//int F2_integrand_A(const int * __restrict  ndim, const PREC* __restrict intv,const int *__restrict ncomp,PREC* __restrict f, void* __restrict p){
int F2_integrand_A(const int  *ndim,const  PREC *intv,const int  *ncomp,PREC * f, void * __restrict p){
	Integrand_kt* param=(Integrand_kt*)p;
	const PREC& beta=intv[0];
	const PREC& kappa_t_prime2=intv[1];
	const PREC& kt2=intv[2];
	PREC val=0;

	val+= (2.0/3.0)*param[0].integrand01( kt2, kappa_t_prime2, beta);
	val+= (4.0/9.0)*param[1].integrand01( kt2, kappa_t_prime2, beta);
	val+= (1.0/9.0)*param[2].integrand02( kt2, kappa_t_prime2, beta);
	val+= (2.0/3.0)*param[0].integrand02( kt2, kappa_t_prime2, beta);
	val+= (4.0/9.0)*param[1].integrand02( kt2, kappa_t_prime2, beta);
	val+= (1.0/9.0)*param[2].integrand02( kt2, kappa_t_prime2, beta);
	//std::cout<<*f<<std::endl;
	if(std::isnan(val)+std::isinf(val)!=0){
		printf("evaluation failure %.5le\n", (double)val);
	}
	*f=val;
	return 0;
}
/////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////
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
	PREC integrand_r(PREC z,PREC r)const{
		PREC val=(Q2)/(2*PI*r) *sigma(r)* psisq_f (z, r);
		return(val);//r^2 comes from photon wave function. just extracted... 2 pi r is angular integration 
	}
};

int F2_integrand_B(const int *__restrict ndim, const PREC *__restrict intv,const int *__restrict ncomp,PREC*__restrict  f, void* __restrict p){
	Integrand_r *param=(Integrand_r*)p;
	PREC z=intv[0];
	PREC r=intv[1];
	
	PREC jac=0;
	change_var(r,jac,R_MIN,R_MAX);
	PREC res=0;
	
	res+=(2.0/3.0)*param[0].integrand_r(z,r);
	res+=(4.0/9.0)*param[1].integrand_r(z,r);
	res+=(1.0/9.0)*param[2].integrand_r(z,r);

	*f=jac*res;
	return(0);
}

int F2_integrand_T(const int * ndim, const PREC intv[],const int * ncomp,PREC f[], void* p){
	*f=intv[0]*intv[1]*intv[1];
	return 0;
}

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

double F2_kt(const PREC x,const  PREC Q2,const  PREC mf2,const PREC (& par)[]){
	static int count;
	std::string type="gbw";
#if R_FORMULA==1
	const int ndim=2;
	Sigma sigma[]={Sigma(type,par) ,Sigma(type,par) ,Sigma(type,par) };
	Integrand_r integrands[]={
		Integrand_r(type,modx(x,Q2,MASS_L2),Q2,MASS_L2,sigma[0]) ,
		Integrand_r(type,modx(x,Q2,MASS_C2),Q2,MASS_C2,sigma[1]) ,
		Integrand_r(type,modx(x,Q2,MASS_B2),Q2,MASS_B2,sigma[2])
       	};
	const int key =13;
#else
	const int ndim=3;
	static Gluon gluon[]={Gluon(type,par) ,Gluon(type,par) ,Gluon(type,par) };
	static Integrand_kt integrands[]={
		Integrand_kt(type, x, Q2, MASS_L2, gluon[0]),
	       	Integrand_kt(type, x, Q2, MASS_C2, gluon[1]),
	       	Integrand_kt(type, x, Q2, MASS_B2, gluon[2])
       	};

	gluon[0].set_par(par);
	gluon[1].set_par(par);
	gluon[2].set_par(par);
	integrands[0].set(x,Q2,MASS_L2);
	integrands[1].set(x,Q2,MASS_C2);
	integrands[2].set(x,Q2,MASS_B2);

	const int key =11;
#endif
	
	const long long int mineval=pow(15,ndim), maxeval=1/pow(INT_PREC/10,2);//use llChure if larger than ~1.0e+9
	const long long int nstart=1.0e+2,nincrease=1.0e+2;
	long long int neval=0;
	//const int mineval=pow(10,ndim), maxeval=1/pow(INT_PREC/10,2);//use llChure if larger than ~1.0e+9
	//const int nstart=1.0e+2,nincrease=1.0e+2;
	//int neval=0;

	const int flag= 0+4*0+8*1+16*0+32*0;
	
	int nregions=0,fail=0;
	//cubareal
	PREC integral[3]={0},error[3]={0},prob[3]={0};
	//int spin=0;
	char statefile[100]="";
 	PREC result=0;
 	//printf("%d: Cuhre x=%.2e Q2=%.2e mf2=%.2e\n",count++, x,Q2,mf2);
	llCuhre(ndim, 1,
	//llTest(ndim, 1,
#if R_FORMULA==1
		&F2_integrand_B,
#else
		&F2_integrand_A,
		//&F2_integrand_A1,
#endif  
		(void*)integrands,
		 1,INT_PREC,INT_PREC/10, flag, mineval,maxeval, key,statefile,NULL, &nregions, &neval,  &fail, integral, error, prob
	);
	//printf("Cuhre-End\n");
	//cubawait(&spin);

 	result=integral[0];
 
 	 if(std::isnan(result)+std::isinf(result)!=0){
 	 	printf("%.3le encountered \n",(double)result);
 	 	return 0;
 	 }
 	// printf("end F2\n");
 	 return ((double)(result));
}



