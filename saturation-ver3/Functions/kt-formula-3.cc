#include<iostream>
#define TEST 1 

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



double change_var(double & var,double &  jac,const double min, const double max){
	return(change_var(var,jac,min,  max,1));
}

class Gluon{
	double sigma_0=0,lambda=0,x_0=0;
	double Q2=0;
	std::string key;
	
	public:
		explicit Gluon(std::string &type,const double (&par)[]){
				//printf("gluon created\n");
				key=type;	
				set_par(par);	
		}
		void set_par(const double(&par)[]){
			if(key=="gbw"){
				//printf("parameters set\n");
				sigma_0 =(double)par[0];
				lambda	=(double)par[1];
				x_0	=(double)par[2];
			}else{
				std::cout<<"unknown model: "<<key<<std::endl;
			}

		}
		~Gluon(){}
		int set_kinem(const double &b){
			//x=a;
			Q2=b;
			return 0;
		}
	public:
		inline double alpha(const double mu2)const{
			return 4.0/(9.0 *log( ((mu2>2*LQCD2)?(mu2):(2.0*LQCD2))/LQCD2));
		}

	public:
		double operator()(const double x,const double k2,double mu2){
			if(x_0<1.0e-5||x_0>1.0e-3){
				return 0;
			}
			if(lambda<0.05||lambda>0.95){
				return 0;
			}
			double Qs2=pow(x_0/x,lambda);
			double val=3.0/(4*PI*PI)*k2/Qs2*exp(-k2/Qs2);
			if(std::isnan(val)==1){
				return(0);
			}
#if ALPHA_RUN==1
			val*=alpha(mu2)/0.15;
#endif
			return (sigma_0*val) ;
		}
};
class Integrand_kt{
	public:
		Integrand_kt(const double a, const double b, const double c, Gluon& af ){
			set_kinem(a,b,c);
			gluptr=&af;
			gluptr->set_kinem(Q2);
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
		~Integrand_kt(){
		
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
			double val=0;//,valtot=0;
			change_var(k2,jac1,0,k2max,1+k2max/pow(Q2,0.25));
			val=integrand(kappa2,k2,beta,phi)+integrand(kappa2,k2,1-beta,phi+PI);
			val*=jac1*jac2*jac3*jac4;

/*			if(kappa2<k2max){
				k2=1-x1*x1;
				change_var(k2,jac1,0,kappa2,1+kappa2/pow(Q2,0.25));
				val=integrand(kappa2,k2,beta,phi)+integrand(kappa2,k2,1-beta,phi+PI);
				val*=jac1*jac2*jac3*jac4;
				valtot+=2*x1*val;
			}else{
				change_var(k2,jac1,0,k2max,1+k2max/pow(Q2,0.25));
				val=integrand(kappa2,k2,beta,phi)+integrand(kappa2,k2,1-beta,phi+PI);
				val*=jac1*jac2*jac3*jac4;
			}
			if(kappa2<k2max){
			
				k2=x1*x1;
				change_var(k2,jac1,kappa2,k2max,1+(k2max-kappa2));
				val=integrand(kappa2,k2,beta,phi)+integrand(kappa2,k2,1-beta,phi+PI);
				val*=jac1*jac2*jac3*jac4;
				valtot+=2*x1*val;
			}
*/
			if(k2>k2max||std::isnan(val)||std::isinf(val)){
				printf("x1= %.3e x2= %.3e x3= %.3e x4= %.3e\n",x1,x2,x3,x4 );
				printf("kappa2=%.3e k2=%.3e beta=%.3e phi=%.3e\n",kappa2,k2,beta,phi );
				printf(" %.3e <beta<  %.3e\n",betamin,betamax );
				printf("%.3e <k2< %.3e\n",0.0,k2max);
				//printf("%.3e <cos[phi]\n",cosphimin );
				printf("%.3e <kappa2\n",kappamax );
				getchar();
			}
			return(val);
		}

	private:
		double mf2=0,Q2=0,x=0;
		double betamin=0,betamax=0,kappamax=0;
		Gluon* gluptr=NULL;

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

		double integrand(double kappa2,double k2,double beta,double phi) const {
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


int F2_integrand_A(const int  *ndim,const  PREC *intv,const int  *ncomp,PREC * f, void * __restrict p){
		Integrand_kt* integrand=(Integrand_kt*)p;
		const PREC beta=intv[0];
		const PREC kappa_t_prime2=intv[1];
		const PREC kt2=intv[2];
		const double phi=intv[3];
		PREC val=0;
		val+= (2.0/3.0)*integrand[0]( kt2, kappa_t_prime2, beta,phi);
		val+= (4.0/9.0)*integrand[1]( kt2, kappa_t_prime2, beta,phi);
		val+= (1.0/9.0)*integrand[2]( kt2, kappa_t_prime2, beta,phi);
		*f=val/(4*PI);
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
	const int ndim=4;
	//create once
	//static F2_integrand F2_integrand_A;
	static Gluon gluon[]={Gluon(type,par) ,Gluon(type,par) ,Gluon(type,par) };
	static Integrand_kt integrands[]={
		Integrand_kt( x, Q2, MASS_L2, gluon[0]),
	       	Integrand_kt( x, Q2, MASS_C2, gluon[1]),
	       	Integrand_kt( x, Q2, MASS_B2, gluon[2])
       	};

	gluon[0].set_par(par);
	gluon[1].set_par(par);
	gluon[2].set_par(par);
	integrands[0].set_kinem(x,Q2,MASS_L2);
	integrands[1].set_kinem(x,Q2,MASS_C2);
	integrands[2].set_kinem(x,Q2,MASS_B2);

	const int key =9;
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

 	result=Q2/(2*PI) *integral[0];
 
 	 if(std::isnan(result)+std::isinf(result)!=0){
 	 	printf("%.3le encountered \n",(double)result);
 	 	return 0;
 	 }
 	// printf("end F2\n");
 	 return ((double)(result));
}














