#include<cmath>
#include<iostream>
#include<fstream>
#include<string>

#include"constants.h"
#include"control.h"
#include"control-default.h"
#include"../Utilities/options.h"
#include"./clenshaw-curtis.hh"


#ifndef ADJOINT
	#define ADJOINT 0
#endif
//double      cyl_bessel_k( double v, double x );
////////////////////////////////////////////////////////////////////////////
//Formulae can be found in J. R. Forshaw, et al. Nucl. Phys. A675, (2000) //
////////////////////////////////////////////////////////////////////////////

extern "C" double SIGMA(double ,double, double,double*,double*);
extern "C" int parameter(double*, double*,double*);
extern "C" int approx_xg(double*);
extern "C" double dclenshaw(double (*)(double*), double,double,double);

double IIM(double r,double x, double Q2){
	double lambda=0.2197;
	double Qs=pow(1.632e-5/x,lambda/2 );
	double alpha=0.85,beta=0.9;
	double gamma=0.7376;
	double kappa=9.9;

	double val=0;
	if(r*Qs<=2){
		val=0.7*pow(r*Qs/2,2*gamma)*exp(-2*pow( log(r*Qs/2),2 )/(kappa*lambda*log(1/x)) );
	}else{
		val=1-exp(-4*alpha*pow( log(beta*r*Qs),2 ) );
	}
	return(val);

}

class impact{
///	private:
//		double z;
	public:
	int index;
	double ep,k,beta,xp,Q2, mf2, kt2,z;
	double *sigpar,*sudpar;
	void set_extern(double a2, double a3,double a4,double a5){
		beta=a2;
		xp=a3;
		Q2=a4;
		mf2=a5;
		}
	void set_z(double a1){
		if(a1>1||a1<0){
			printf("set_z:: can not set  z = %f\n",z);
		}
		z=a1;
		ep=sqrt(fabs(z*(1-z)*Q2+mf2));
		k=sqrt(fabs(z*(1-z)*Q2*(1-beta)/beta -mf2));
		//std::cout<<"z= "<<z<<"\tk= "<<k<<"\tep= "<<ep<<std::endl;
	}
	void current(){
		std::cout<<"beta: "<<beta<<"\t xp: "<<xp<<"\t Q2: "<<Q2<<"\t mf2: "<<mf2<<std::endl;
		std::cout<<"z: "<<z<<"\t "<<"k: "<<k<<"\t "<<"ep: "<<ep<<std::endl;
	}
};
/////////////////////////////////////
//////  stores global variables /////
static struct impact diff_param, diff_param_beta;
//Clenshaw_Curtis integral(32);
////////////////////////////////////

double phi_integrand(double *R){
	double r=*R;//change of variable
	double jac=pow(1-r,-2);
	r=r/(1-r);

	double k=diff_param.k, ep=diff_param.ep,  index=diff_param.index;
	//std::cout<<"z= "<<diff_param.z<<"\t"<<"kr= "<<k<<" * "<<r<<"\t"<<"ep r= "<<ep<<" * "<<r<<std::endl;
	double x=diff_param.xp;//Q2*(1/xp-1);
	

	double val=r*std::cyl_bessel_j(index,k*r)*std::cyl_bessel_k(index, ep*r);//	*(diff_param.func)(r,x,Q2,sigpar);
	if(isnan(r)==1){
		printf("%f passed from integral func, makes r=%f\n",*R,r);
	}
	//val*=SIGMA(r,x,diff_param.Q2,diff_param.sigpar,diff_param.sudpar);
	val*=diff_param.sigpar[0]*IIM(r,x,diff_param.Q2);
	val*=jac;

	if((isnan(val)+isinf(val))!=0){
		diff_param.current();
		printf("val=%.3e,r= %.3e\n",val,*R);
	}
	return(val);
}

Clenshaw_Curtis phi_integrator(16);

double phi(int index,double z){
	diff_param.index=index;
	diff_param.set_z(z);

	//double val=phi_integrator.integrate(&phi_integrand,1.0e-5,0.99,1.0e-6);
	double val=dclenshaw(&phi_integrand,1.0e-5,0.99,1.0e-6);
	return(val);
}

double phi2_integrand_u(double *U){
	double u=*U;//change of variable
	double jac=pow(1-u,-2);
	u=u/(1-u);
	
	double kt=std::sqrt(diff_param.kt2) , z=diff_param.z, xp=diff_param.xp, Q2=diff_param.Q2;
	int index=2;
	double x=xp;//Q2*(1/xp-1);
	
	double val=u*std::cyl_bessel_k(index,std::sqrt(z/(1-z))*u)*std::cyl_bessel_j(index, u );//	*(diff_param.func)(r,x,Q2,sigpar);
	if(isnan(u/kt)==1){
		printf("%f passed from integral func, makes u=%f\t kt2=%f\n ",*U,u,diff_param.kt2);
	}
	//double sigma=SIGMA(u/kt,x,Q2,diff_param.sigpar,diff_param.sudpar);
	double sigma=diff_param.sigpar[0]*IIM(u/kt,x,Q2);
#if ADJOINT==0
	val*=9.0/4.0*sigma;
#elif ADJOINT==1
	sigma/=( diff_param.sigpar[0]);
	//val*=(1-pow(1-sigma, 9.0/4.0))*diff_param.sigpar[0];
	val*=(2*sigma-sigma*sigma)*( diff_param.sigpar[0]);
#endif
	val*=jac;
	return(val);
}

Clenshaw_Curtis phi2u_integrator(16);
double phi2_integrand_kt(double *K){
	double kt2=pow(*K,2);
	double jac=2**K;
	
	if(kt2<0){
		diff_param.current();
		getchar();
		//return(0);
	}
	diff_param.kt2=kt2;

	//double phi2=phi2u_integrator.integrate(&phi2_integrand_u,1.0e-5,0.98, 1.0e-4);
	double phi2=dclenshaw(&phi2_integrand_u,1.0e-5,0.99,1.0e-5);

	double z=diff_param.z, Q2=diff_param.Q2;
	double beta=diff_param.beta;
	double val=phi2*phi2*std::log((1-z)*Q2/(kt2)/**(z-beta)/beta*/);
	return(jac*val);
}

////////////////////////////////////////////////////////////////////
double FD_L_integrand(double *Z){
	//std::cout<<"FDL z= "<<*Z<<std::endl;
	double z=*Z;
	//double z=pow(*Z,2)/2;
	double jac=1;//*Z;
	if(z<0||z>1){
		printf("FD_L_integrand:: z = %.5e\n",z);
	}
	double phi0=phi(0,z);
	double val=phi0*phi0*pow(z * (1-z),3);
	//if(isinf(val)==1||isnan(val)==1){
	//	printf("FD_L_integrand:: z=%.3e phi0=%.3e val=%.3e\n",*Z,phi0,val);
	//}
	return(jac*val);
}
double FD_T_integrand(double *Z){
	//std::cout<<"FDT z= "<<*Z<<std::endl;
	double z=*Z;
	//double z=pow(*Z,2)/2;
	double jac=1;//*Z;
	if(z<0||z>1){
		printf("FD_T_integrand:: z = %.5e\n",z);
	}
	double phi0=phi(0,z);
	double phi1=phi(1,z);
	//!!! NOTE !!!, in the phi function diff_parameters k and ep are modified as they depend on z.
	double val=pow(diff_param.ep,2)*(z*z+(1-z)*(1-z))*phi1*phi1;
	val+=diff_param.mf2*phi0*phi0;
	val*=z*(1-z);
	//if(isinf(val)==1||isnan(val)==1){
	//	printf("FD_T_integrand::  z=%.3e phi0=%.3e, phi1=%.3e val=%.3e\n",*Z,phi0,phi1,val);
	//}
	return(jac*val);
}	
Clenshaw_Curtis phi2_integrator(16);
double FD_g_integrand(double *Z){
	double z=*Z;
	if(z<0||z>1){
		printf("FD_g_integrand:: z = %.5e\n",z);
	}
	//diff_param.index=index;
	diff_param.set_z(z);
	double beta=diff_param.beta;
	double val,min=1.0e-7,max=(1-z)*diff_param.Q2;//*(z-beta)/beta;
	max=sqrt(max);

	if(min>max){
		return(0);
	}

	if(min<0){
		printf("%f !!!!!\n",min);
		getchar();
	}
	//val=phi2_integrator.integrate(&phi2_integrand_kt, min,max,1.0e-4);
	
	val=dclenshaw(&phi2_integrand_kt,min,max,1.0e-4);
	val*=(pow(1-beta/z,2)+pow(beta/z,2) )/pow(1-z,3);
	if(isinf(val)==1||isnan(val)==1){
		printf("FD_g_intregrand:: z=%.3e val=%.3e\n",*Z,val);
	}
	return(val);
}




///////////////////////////////////////////////////////////////////////////////////////////////
// Small beta
///////////////////////////////////////////////////////////////////////////////////////////////
Clenshaw_Curtis AIntegrator(16), AIntegrator2(16);
struct A_param{
	double xp;
	double v,r;	
}A_var;

double S(double r,double xp){
	//double 	val=1 - SIGMA(r,xp,diff_param_beta.Q2,diff_param_beta.sigpar,diff_param_beta.sudpar)/( diff_param_beta.sigpar[0] );
	double 	val=1 - IIM(r,xp,diff_param_beta.Q2);
	return(val);
}
double A_integrand(double u,double v, double r,double xp){
	double val=0;
	val+=u*u*pow( 
			S( r/u ,xp )*S( (1.0/u - 1+2*v)*r,xp )-S(r,xp) 
		,2);
	val+=pow( 
			S( r*u ,xp )*S( (1.0 - u+2*u*v)*r,xp )-S(r,xp) 
		,2);
	val/=( u*(1-u+2*u*v)*sqrt(v*(1-v)*(1+u*v)*(1-u+u*v)) );
	val*=2;
	return(val);
}

double AIntegrand1(double *U){
	double val=A_integrand(*U,A_var.v,A_var.r, A_var.xp);
	return(val);
}
double AIntegrand2(double *V){
	A_var.v=*V;
	//integrate over u;
	double val=AIntegrator.integrate(&AIntegrand1,1.0e-8,1-1.0e-8,1.0e-3);
	return(val);
}
double A(double r,double xp){
	A_var.xp=xp;
	A_var.r=r;
	double val=AIntegrator2.integrate(&AIntegrand2,1.0e-8,1-1.0e-8,1.0e-2);
	return(val);
}


Clenshaw_Curtis xF_Integrator1(16), xF_Integrator2(16);
struct phi_param{
	double r;
}phi_T_param;

double phi_T(double z, double r, double mf2){
	double val=0;
	double ep=diff_param_beta.ep;
	val+=(z*z+(1-z)*(1-z))*ep*ep*pow(std::cyl_bessel_k(1,ep*r),2);
	val+=mf2*pow(std::cyl_bessel_k(0,ep*r),2);
	return(val);
}
double xF_beta_Integrand1(double *Z){
	diff_param_beta.set_z(*Z);
	double val=phi_T(*Z,phi_T_param.r,diff_param_beta.mf2);
	//double val=xF_Integrator1.integrate(&xF_beta_Integrand1, 1.0e-8,100,1.0e-3);
	return(val);
}

double xF_beta_Integrand2(double *R){
	phi_T_param.r=*R;
	double val=xF_Integrator1.integrate(&xF_beta_Integrand1, 1.0e-8,1-1.0e-8,1.0e-3);
	val*=*R* A(*R,diff_param_beta.xp);
	return(val);

}
//double xF_beta_Integrand2(double *R){
//	double val=*R* A(*R,diff_param_beta.xp)*phi_T(diff_param_beta.z,*R,diff_param_beta.mf2);
//	return(val);
//}
	

double xF_beta(double xp,double Q2,double mf2){
	getchar();
	//diff_param.set_extern(0,xp,Q2,mf2);
	double val=xF_Integrator2.integrate(&xF_beta_Integrand2, 1.0e-8,100,1.0e-2);
	val*=Q2*3.0/(4.0*PI);

	return(val);
}

//////////////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////////////
Clenshaw_Curtis xF_beta0_Integrator1(16), xF_beta0_Integrator2(16);
double xF_beta0_Integrand1(double *R){
	double r=*R;
	double xp=diff_param.xp;
	double k=sqrt(diff_param.kt2);
	//double N=SIGMA(r,xp,diff_param.Q2,diff_param.sigpar,diff_param.sudpar)/( diff_param.sigpar[0] );
	double N=IIM(r,xp,diff_param.Q2);
	double val=std::cyl_bessel_j(2,r*k)*(2*N-N*N)/r;
	return(val);
}
double xF_beta0_Integrand2(double *K){
	double k=*K;
	diff_param.kt2=k*k;
	double Q2=diff_param.Q2;

	double val=xF_beta0_Integrator1.integrate(&xF_beta0_Integrand1,1.0e-8,100,1.0e-3);
	val=val*val;
	val*=2*k*log(Q2/(k*k));
	return(val);
}
double xF_beta0(double xp,double Q2,double mf2){
	//diff_param.set_extern(0,xp,Q2,mf2);
	double val=xF_beta0_Integrator2.integrate(&xF_beta0_Integrand2,1.0e-8,sqrt(Q2),1.0e-3);
	return(val);
}

//////////////////////////////////////////////////////////////////////////////////////////
//
//
Clenshaw_Curtis FD3_integrator(16);
double xFD_LT(int pol){
	//!!!!!!!!  Run diff_param.set_extern(beta,xp,Q2,mf2) before use !!!!!!!!!!!
	const double hc22=pow(0.3893,2), BD=5.591, alpha_s=0.25;//hc2, GeV to mb conversion.
	double factor, min,max , beta_factor;
	double (*funcptr)(double*);
	getchar();
	switch(pol){
		case 't':
			factor=2*3*pow(diff_param.Q2,2)/(hc22*128*pow(PI,4)*diff_param.beta*BD);
			funcptr=&FD_T_integrand;
			max=0.5;
			//Half the regon. cf. GBW int over k for Q2(1-beta)/(4*beta)>k2>0.
			// sym. z<->1-z. so factor 2.
			min=(1-std::sqrt(1-4*(diff_param.mf2*diff_param.beta/(diff_param.Q2*(1-diff_param.beta))) ) )/2;
			//max=sqrt(2*max);
			min=((min>1.0e-10)?(min):(1.0e-10));
			break;

		case 'l':
			factor=2*3*pow(diff_param.Q2,3)/(hc22*32*pow(PI,4)*diff_param.beta*BD);
			funcptr=&FD_L_integrand;
			max=0.5;
			min=(1-std::sqrt(1-4*(diff_param.mf2*diff_param.beta/(diff_param.Q2*(1-diff_param.beta))) ) )/2;
			//max=sqrt(2*max);
			//min=sqrt(2*min);
			min=((min>1.0e-10)?(min):(1.0e-10));
			break;
		case 'g':
#if ADJOINT==0
			factor=81*diff_param.beta*alpha_s/(512*pow(PI,5)*hc22*BD);
#elif ADJOINT==1	
			//printf("beta_factor\n");
			beta_factor=xF_beta(diff_param.xp,diff_param.Q2,diff_param.mf2)/ xF_beta0(diff_param.xp,diff_param.Q2,diff_param.mf2);
			//beta_factor=1;
			//printf("FACTOR\t %.3e\n",beta_factor);

			factor=diff_param.beta*alpha_s/(32*pow(PI,5)*hc22*BD);
			factor*=beta_factor;
#endif
			funcptr=&FD_g_integrand;
			max=1-1.0e-8;
			min=diff_param.beta+1.0e-8;
			break;

		default:
			printf("error:: choose polarizaion\n");
		}
	//std::cout<<"z integral ["<<min<<", "<<max<<"]"<<std::endl;
	double res=FD3_integrator.integrate(funcptr,min,max,1.0e-2);
	//res*=2;//see the comm. above.
	std::cout<<"pol="<<(char)pol<<"\tIntegral= "<<res;

	res*=factor;
	res*=2.0/3.0;//sum of charges uds;
	std::cout<<"\tResult=  "<<res<<std::endl;
	return(res);
}
/////////////////////////////////////////////////////////////////////////////////////////
double xFD_T(double beta,double xp,double Q2,double mf2){
	diff_param.set_extern(beta,xp,Q2,mf2);
	double res=xFD_LT('t');
	return(res);
}
double xFD_L(double beta,double xp,double Q2,double mf2){
	diff_param.set_extern(beta,xp,Q2,mf2);
	double res=xFD_LT('l');
	return(res);
}
double xFD_g(double beta,double xp,double Q2,double mf2){
	diff_param.set_extern(beta,xp,Q2,mf2);
	double res=xFD_LT('g');
	return(res);
}
///////////////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////////////


int main(int argc,char** argv){
	char input_file_name[250];
	strcpy(input_file_name,argv[1]);
	char output_file_name[250];
	strcpy(output_file_name,argv[2]);
	double Q2=atof(argv[3]);
	double beta=atof(argv[4]);
	double xmin=atof(argv[5]),xmax=atof(argv[6]);
	printf("%.3e %.3e %.3e %.3e\n",Q2,beta,xmin,xmax); 	

	//read_options(argc,argv,&OPTIONS);
	//FILE* infile=fopen(OPTIONS.input_file_name,"r");
	FILE* infile=fopen(input_file_name,"r");
	double param_arr[20],sigpar[10],sudpar[10];
	read_parameters(infile,param_arr);
	parameter(param_arr,sigpar,sudpar);

	sigpar[0]=27.36;
	diff_param.sigpar=sigpar;
	diff_param.sudpar=sudpar;
	diff_param_beta.sigpar=sigpar;
	diff_param_beta.sudpar=sudpar;

       	fclose(infile);
	//FILE* controlfile=fopen(argv[argc-1],"r" );
	//double Q2,beta,xmin,xmax;
	//fscanf(controlfile,"%lf\t%lf\t%lf%lf",&beta,&Q2,&xmin,&xmax );
	//fclose(controlfile);
	//Clenshaw_Curtis integral(16);

	//phi2u_integrator.Print=1;
	phi2u_integrator.DIV=4;
#if ((MODEL==1)||(MODEL==3))
	approx_xg(sigpar+1);
#endif

	//beta=OPTIONS.beta;
	//Q2=OPTIONS.Q2;
	//double mf2=MASS_L2;
	double mf2=0.0196;
	beta*=(1+4*mf2/Q2);
	double xp;
	double val;
	std::fstream file;
	
	//char name[150];
	//sprintf(name,"%s-%.4f-%.3f.txt",OPTIONS.output_file_name,beta,Q2);

	//file.open(OPTIONS.output_file_name,std::fstream::out);
	file.open(output_file_name,std::fstream::out);
	printf("%s\n", output_file_name);
	char type[3]={'t','l','g'};
	//xmin=(OPTIONS.xmin )/2;
	//xmax=(OPTIONS.xmax )*2;
	xmin=(xmin )/2;
	xmax=(xmax )*2;
	std::cout<<"Q2= "<<Q2<<"\tbeta= "<<beta<<std::endl;
	double y;

	for(int i=0;i<5;i++){
		xp=pow(10,log10(xmin)  + log10(xmax/xmin )*((double)i)/(5-1))/beta;
		printf("%.5e\t %.5e\t %.5e\n",xp,beta,Q2);
		diff_param.set_extern(beta,xp,Q2,mf2);
		diff_param_beta.set_extern(beta,xp,Q2,mf2);
		val=0;
		for(int j=0;j<3;j++){
			if(type[j]=='l'){
				y=Q2/(318.12*318.12*xp*beta);//H1 sqrt(s)=318.12 GeV
				val+=2*((1-y)/(1+(1-y)*(1-y)))*xFD_LT(type[j]);

			}else{
				val+=xFD_LT(type[j]);
			}
		}
		std::cout<<"xp= "<<xp<<"\tval= "<<val<<std::endl;
		file<<xp<<"\t"<<val<<std::endl;
	}
	file.close();

	return(0);
}
