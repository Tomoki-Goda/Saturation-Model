#include<iostream>

#include<cuba.h>
#include<cmath>
#include"./control.h"
#include"./control-default.h"
#include"./constants.h"




extern PREC INT_PREC;
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
	if(std::isnan(x+Q2+mf2+mf2/Q2)+std::isinf(x+Q2+mf2+mf2/Q2)!=0){
				printf("modx:: Q2=%.3le x=%.3le mf2=%.3le\n",(double)Q2,(double)x,(double)mf2);
				getchar();
	}
#if MODX==1
	return( (x*(1+4*mf2/Q2)));
#else 
	return( x);
#endif
}

int change_var(PREC* var,PREC* jac,const PREC min, const PREC max){
	PREC scale=(max-min);
	PREC a1=*var;
	PREC a2=1-a1;
	PREC den=scale*a2+a1;
	*jac=pow(scale/den,2);
	*var=(min*(scale)*a2+max*a1)/den;
	return 0;
}

class Gluon{
	PREC sigma_0,lambda,x_0;
	PREC Q2;
	std::string key;
	
	public:
		Gluon(std::string type,const double*par){
			if(type=="gbw"){
					sigma_0 =(PREC)par[0];
					lambda	=(PREC)par[1];
					x_0	=(PREC)par[2];
					key=type;
			}else{
					std::cout<<"unknown model: "<<type<<std::endl;
			}
		
		}
		~Gluon(){
		}
		int set_kinem(const PREC b){
			//x=a;
			Q2=b;
			return 0;
		}
	private:
		inline PREC alpha(const PREC mu2){
			//const PREC b0=((11.0*3.0-3.0*2.0)/12.0);
			//const PREC val=1.0/(b0 *log( ((mu2>2*LQCD2)?(mu2):(2*LQCD2))/LQCD2));
			return 4.0/(9.0 *log( ((mu2>2*LQCD2)?(mu2):(2*LQCD2))/LQCD2));
		}

	public:
		PREC aF(const PREC x,const PREC k2){
			//PREC lambda=par[1];
			//PREC x0=par[2];

			PREC Qs2=pow(x_0/x,lambda);
			PREC val=3.0/(4*PI*PI)*k2/Qs2*exp(-k2/Qs2);



			if(std::isnan(val)==1){
				return(0);
			}
			return (sigma_0*val) ;
		}
};

class Integrand_kt{
		Gluon * gluptr;
		PREC x,Q2,mf2;
		
	public:
		Integrand_kt(const std::string type,const PREC a,const PREC b,const PREC c,Gluon* gluon){
			gluptr=gluon;
			set_kinem(a,b,c);
			gluon->set_kinem(b);
		}
		~Integrand_kt(){
		}
		
		int set_kinem(const PREC a,const PREC b,const PREC c){
			x=a;
			Q2=b;
			mf2=c;
			return 0;
		}
		
	private:
		int I_array(const PREC beta,const  PREC kappa_t_prime2,const PREC kt2, PREC *I)const{
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



		//static 
		inline PREC inv_z(const PREC beta, const PREC kappa_t_prime2,const PREC kt2)const{
			return( 1+(kappa_t_prime2+mf2)/(beta*(1-beta)*Q2)+kt2/Q2);
		}

		

	public:
		PREC  integrand(PREC kt2, PREC kappa_t_prime2, PREC beta)const{
			if((1-x)/x *Q2-4*mf2<=0){
				return(0);
			}
			PREC jac1,jac2,jac3;
			
			const PREC ktmax=(1-x)/x *Q2-4*mf2;
			change_var(&kt2, &jac1,0,ktmax);
			//PREC kprimemax=(1-x)/x *Q2/4-mf2;
			const PREC kprimemax=ktmax/4 ;//just because it is...
			change_var(&kappa_t_prime2, &jac2,0,kprimemax);
			PREC betamin=sqrt(1-4*(mf2/((1-x)/x*Q2) )), betamax ;
			betamax=(1+betamin)/2;
			betamin=(1-betamin)/2;
			change_var(&beta, &jac3,betamin,betamax);
			if(std::isnan(jac1*jac2*jac3)==1){
				printf("Q2= %.3le, x=%.3le, mf2=%.3le sqrt : %.3le jac1 %.3le jac2 %.3le jac3 %.3le  \n",(double)Q2,(double)x,(double)mf2,(double)(1-4*(mf2/((1-x)/x*Q2) )) ,(double)jac1,(double) jac2,(double)jac3);
			}
			const PREC xz=x*inv_z(beta,kappa_t_prime2,kt2) ;
			if(1-xz<0){
				return 0;		
			}
			
			PREC val=0;
			PREC I[4];
			I_array(beta, kappa_t_prime2, kt2, I);

			val+=(beta*beta+pow(1-beta,2))*(I[0]-I[1]);
			val+=(mf2+4*Q2*beta*beta*pow(1-beta,2) )*(I[2]-I[3]);
			val*=gluptr->aF(xz, kt2 );
			
#if ALPHA_RUN==1
			val*=alpha(kappa_t_prime2+kt2+mf2+1)/0.2;
			//printf("%.5le\n", val);
#endif
			return(Q2/(2*PI) *jac1*jac2*jac3* val/kt2 );	
			//return( jac1*jac2*jac3* F2_integrand(beta,kappa_t_prime2,kt2,Q2,mf2,x,par )/kt2 );
		}
};


int F2_integrand_A(const int *ndim, const PREC* intv,const int *ncomp,PREC* f, void* p){
	Integrand_kt* param=(Integrand_kt*)p;
	PREC beta=intv[0];
	PREC kappa_t_prime2=intv[1];
	PREC kt2=intv[2];
	PREC val=0;

	val+= (2.0/3.0)*param[0].integrand( kt2, kappa_t_prime2, beta);
	val+= (4.0/9.0)*param[1].integrand( kt2, kappa_t_prime2, beta);
	val+= (1.0/9.0)*param[2].integrand( kt2, kappa_t_prime2, beta);
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

		Sigma(std::string  type , const PREC* par){//maybe use struct pointer for parameters.
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
		Integrand_r(std::string model,PREC x,PREC Q2,PREC mf2,Sigma* sig){
			sigma_ptr=sig;
			set_kinem(x,Q2,mf2);
			//printf("int created %.3e\n",(double)integrand_r(0.1,0.1));
		};
	private:
		PREC x, Q2, mf2;
		int set_kinem(PREC a,PREC b,PREC c){
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

int F2_integrand_B(const int *ndim, const PREC *intv,const int *ncomp,PREC* f, void* p){
	//PREC **param=(PREC**)p;
	
	Integrand_r *param=(Integrand_r*)p;
	PREC z=intv[0];
	PREC r=intv[1];
	
	PREC jac=0;
	change_var(&r,&jac,R_MIN,R_MAX);
	PREC res=0;
	
	res+=(2.0/3.0)*param[0].integrand_r(z,r);
	res+=(4.0/9.0)*param[1].integrand_r(z,r);
	res+=(1.0/9.0)*param[2].integrand_r(z,r);

	*f=jac*res;
	return(0);
}

int F2_integrand_T(const int *ndim, const PREC intv[],const int *ncomp,PREC f[], void* p){
	*f=intv[0]*intv[1]*intv[1];
	return 0;
}

double F2_kt(const PREC x,const  PREC Q2,const  PREC mf2,const PREC* par){
	//printf("F2\n");
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
	
	const long long int mineval=pow(10,ndim), maxeval=1/pow(INT_PREC/10,2);//use llChure if larger than ~1.0e+9
	const long long int nstart=1.0e+2,nincrease=1.0e+2;

	const int flag= 0+4*0+8*1+16*0+32*0;
	long long int neval=0;
	int nregions=0,fail=0;
	//cubareal
	PREC integral[3]={0},error[3]={0},prob[3]={0};
	
	char statefile[100]="";
 	PREC result;
 	
 	//printf("integ\n");
 	//getchar();
 	
 	//PREC intv[]={0.5, 0.5};
 	//int ncomp=0;
 	//PREC f=0;
 	//F2_integrand_B(&ndim,intv,&ncomp,&f, (void*)(integrands));
 	//printf( "%.3e\n", (double)f);
 	//printf( "%.3e\n", (double)integral[0]);
 	//Cuhre(ndim, 1,&F2_integrand_T,(void*)error, 1,1,1, 1, 10,100, 9,"",NULL, &nregions, &neval,  &fail, integral, error, prob);
 	//printf( "%.3e\n", (double)integral[0]);
	llCuhre(ndim, 1,
#if R_FORMULA==1
		&F2_integrand_B,
#else
		&F2_integrand_A,
#endif  
		(void*)integrands,
		 1,INT_PREC,INT_PREC/10, flag, mineval,maxeval, key,statefile,NULL, &nregions, &neval,  &fail, integral, error, prob
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
 	 	printf("%.3le encountered \n",(double)result);
 	 	return 0;
 	 }
 	// printf("end F2\n");
 	 return ((double)result);
}



