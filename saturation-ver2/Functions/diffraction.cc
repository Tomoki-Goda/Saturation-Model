#include<cmath>
#include<iostream>
#include<fstream>

#include"constants.h"
#include"control.h"
#include"control-default.h"
#include"../Utilities/options.h"


//double      cyl_bessel_k( double v, double x );
////////////////////////////////////////////////////////////////////////////
//Formulae can be found in J. R. Forshaw, et al. Nucl. Phys. A675, (2000) //
////////////////////////////////////////////////////////////////////////////

extern "C" double SIGMA(double ,double, double,double*,double*);
extern "C" int parameter(double*, double*,double*);
extern "C" int approx_xg(double*);
extern "C" double dgquad_(double (*func)(double*),double*, double *, int*);
extern "C" double dgauss_(double (*)(double*),double*,double*,double*);
extern "C" int dadapt_(double (*)(double *), double*,double*,int*,double*,double*,double*,double*);


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
		//std::cout<<beta<<"\t "<<xp<<"\t "<<Q2<<"\t "<<mf2<<std::endl;
		}
	void set_z(double a1){
		z=a1;
		ep=sqrt(fabs(z*(1-z)*Q2+mf2));
		k=sqrt(fabs(z*(1-z)*Q2*(1-beta)/beta -mf2));
		//std::cout<<"z= "<<z<<"\tk= "<<k<<"\tep= "<<ep<<std::endl;
	}
};
/////////////////////////////////////
//////  stores global variables /////
static struct impact diff_param;
////////////////////////////////////
#define INTSTR 1
double integrate(double (*func)(double*), double min,double max ){
	/////////////quad//////////////
#if INTSTR==0
	int N=96;
	double val=dgquad_(func,&min,&max,&N);
#elif INTSTR==1
	int N=96;
	int count=30;
	double high=min,low=min,step=(max-min)/count;
	double val=0;
	for (int i=0;i<count;i++){
		high+=step;
		val+=dgquad_(func,&low,&high,&N);
		low=high;
	}
#elif INTSTR==2
	double rel=1.0e-5; abs=0;
	int seg=1;
	int val,err;
	dadapt_(func, &min,&max,&seg,&rel,&abs,&val,&err);
#elif INTSTR==3
	double eps=1.0e-5;
	double val=dgauss_(func,&min,&max,&eps);
#endif
	return(val);	
}

double phi_integrand(double *R){
	double r=*R;//change of variable
	double jac=pow(1-r,-2);
	r=r/(1-r);

	double k=diff_param.k, ep=diff_param.ep,  index=diff_param.index;
	//std::cout<<"z= "<<diff_param.z<<"\t"<<"kr= "<<k<<" * "<<r<<"\t"<<"ep r= "<<ep<<" * "<<r<<std::endl;
	double x=diff_param.xp;//Q2*(1/xp-1);
	
	double val=r*std::cyl_bessel_j(index,k*r)*std::cyl_bessel_k(index, ep*r);//	*(diff_param.func)(r,x,Q2,sigpar);
	val*=SIGMA(r,x,diff_param.Q2,diff_param.sigpar,diff_param.sudpar);
	val*=jac;
	return(val);
}

double phi(int index,double z){
	diff_param.index=index;
	diff_param.set_z(z);

	//double val=0,min=1.0e-5,max=0.99;
	//double eps=1.0e-6;
	//val=dgauss_(&phi_integrand,&min,&max,&eps);
	double val=integrate(&phi_integrand,1.0e-5,0.99);
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
	val*=SIGMA(u/kt,x,Q2,diff_param.sigpar,diff_param.sudpar);
	val*=jac;
	return(val);
}

double phi2_integrand_kt(double *K){
	double kt2=*K;
	diff_param.kt2=kt2;
	//double min=1.0e-5, max=0.99;
	//int N=96;
	//double phi2=dgquad_(&phi2_integrand_u,&min,&max,&N );
	//double phi2=0;
	//double rel=1.0e-5, abs=0,err;
	//int seg=1;
	//dadapt_(&phi2_integrand_u,&min,&max,&seg,&rel,&abs,&phi2,&err);
	//double eps=1.0e-6;
	double phi2=integrate(&phi2_integrand_u,1.0e-5,0.99);
	double z=diff_param.z, Q2=diff_param.Q2;
	double val=phi2*phi2*std::log((1-z)*Q2/(kt2));
	return(val);
}

////////////////////////////////////////////////////////////////////
double FD_L_integrand(double *Z){
	//std::cout<<"FDL z= "<<*Z<<std::endl;
	double z=*Z;
	double phi0=phi(0,z);
	double val=phi0*phi0*pow(z * (1-z),3);
	return(val);
}
double FD_T_integrand(double *Z){
	//std::cout<<"FDT z= "<<*Z<<std::endl;
	double z=*Z;
	double phi0=phi(0,z);
	double phi1=phi(1,z);
	//!!! NOTE !!!, in the phi function diff_parameters k and ep are modified as they depend on z.
	double val=pow(diff_param.ep,2)*(z*z+(1-z)*(1-z))*phi1*phi1;
	val+=diff_param.mf2*phi0*phi0;
	val*=z*(1-z);
	return(val);
}	
double FD_g_integrand(double *Z){
	double z=*Z;
	//diff_param.index=index;
	diff_param.set_z(z);
	double beta=diff_param.beta;
	int N=96;
	double val,min=1.0e-4,max=(1-z)*diff_param.Q2;
	val=integrate(&phi2_integrand_kt, min,max);
	//double rel=1.0e-5, abs=0,err;
	//int seg=1;
	//dadapt_(&phi2_integrand_kt,&min,&max,&seg,&rel,&abs,&val,&err);

	//double eps=1.0e-5;
	//val=dgauss_(&phi2_integrand_kt,&min,&max,&eps);
	
	val*=(pow(1-beta/z,2)+pow(beta/z,2) )/pow(1-z,3);
	return(val);
}

////////////////////////////////////////////////////////////////////
double xFD_LT(int pol){
	//!!!!!!!!  Run diff_param.set_extern(beta,xp,Q2,mf2) before use !!!!!!!!!!!
	const double hc22=pow(0.3893,2), BD=6.0, alpha_s=0.2;//Unknown hc2 inherited.
	double factor, min,max;
	double (*funcptr)(double*);
	switch(pol){
		case 't':
			factor=2*3*pow(diff_param.Q2,2)/(hc22*128*pow(PI,4)*diff_param.beta*BD);
			funcptr=&FD_T_integrand;
			max=0.5;
			//Half the regon. cf. GBW int over k for Q2(1-beta)/(4*beta)>k2>0.
			// sym. z<->1-z. so factor 2.
			min=(1-std::sqrt(1-4*(diff_param.mf2*diff_param.beta/(diff_param.Q2*(1-diff_param.beta))) ) )/2;
			break;

		case 'l':
			factor=2*3*pow(diff_param.Q2,3)/(hc22*32*pow(PI,4)*diff_param.beta*BD);
			funcptr=&FD_L_integrand;
			max=0.5;
			min=(1-std::sqrt(1-4*(diff_param.mf2*diff_param.beta/(diff_param.Q2*(1-diff_param.beta))) ) )/2;

			break;
		case 'g':
			factor=81*diff_param.beta*alpha_s/(512*pow(PI,5)*hc22*BD);
			funcptr=&FD_g_integrand;
			max=1;
			min=diff_param.beta;
			break;

		default:
			printf("error:: choose polarizaion\n");
		}

	//int N=96;
	//double res=dgquad_(funcptr,&min,&max, &N);
	//std::cout<<"min= "<<min<<"\tmax= "<<max<<std::endl;
	//double res=0;
	/*double high=min,low=min,step=(max-min)/10;
	for(int i=0;i<10;i++){
		high=low+step;
		res+=dgquad_(funcptr,&low,&high,&N);
		low=high;
	}*/
	//double err =0, abs=1.0e-5,rel=1.0e-5;
	//int seg=3;
	//dadapt_(funcptr,&min,&max,&seg,&rel,&abs,&res,&err);

	double eps=1.0e-3;
	double res=dgauss_(funcptr,&min,&max,&eps);
	//res*=2;//see the comm. above.
	std::cout<<"pol="<<pol<<"\t"<<res;

	res*=factor;
	res*=2.0/3.0;//sum of charges uds;
	std::cout<<"  "<<res<<std::endl;
	return(res);
}

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



int main(int argc,char** argv){
	read_options(argc,argv,&OPTIONS);
	FILE* infile=fopen(OPTIONS.input_file_name,"r");
	double param_arr[20],sigpar[10],sudpar[10];
	read_parameters(infile,param_arr);
	parameter(param_arr,sigpar,sudpar);
	diff_param.sigpar=sigpar;
	diff_param.sudpar=sudpar;

       	fclose(infile);
	
#if ((MODEL==1)||(MODEL==3))
	approx_xg(sigpar+1);
#endif

	double beta=0.4;
	double Q2=8.5;
	double mf2=MASS_L2;
	double xp;
	double val;
	std::fstream file;
	file.open(OPTIONS.output_file_name,std::fstream::out);
	char type[3]={'t','l','g'};

	for(int i=0;i<50;i++){
		xp=pow(10,-4+3*((double)i)/50);
		diff_param.set_extern(beta,xp,Q2,mf2);
		val=0;
		for(int j=0;j<3;j++){
			val+=xFD_LT(type[j]);
		}
		std::cout<<xp<<"\t"<<val<<std::endl;
		file<<xp<<"\t"<<val<<std::endl;
	}
	file.close();

	return(0);
}
