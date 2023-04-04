#include<iostream>
#include<fstream>
#include<cmath>

#include<fstream>
#include<vector>
#include<chrono>
#include<cuba.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include<gsl/gsl_dht.h>
#include<gsl/gsl_deriv.h>
#include<gsl/gsl_chebyshev.h>
//#include"Functions/r-formula.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sum.h>
#define N_CHEB 20
#include"Functions/gluons.hh"
#include"Functions/Kahn.hh"
#include"Functions/clenshaw.hh"
#include"Functions/control-default.h"
#include"Functions/constants.h"
#include"Functions/Levin.hh"


///////////////////////////////////////////////////////
//CUBA
/////////////////////////////////////////////////////
/*
class testfunc_cuba{
		double k=0;
		double x=0;
		double rmin=0,rmax=0;
	
	public:
		int set_var(double k,double x, double rmin, double rmax){
			this->k=k;
			this->x=x;
			this->rmin=rmin;
			this->rmax=rmax;
			return 0;
		}
		double operator()(double r){
			r=rmin+rmax*r;

			double qs2=pow(2.0e-4/x,0.2); 
			
			return(std::cyl_bessel_j(0,r*k )* (1-exp(-r*r *qs2/4))/r);
		}

};

template <typename T>int integrand_cuba(const int *__restrict ndim, const double  *__restrict intv,const int *__restrict ncomp,double *__restrict  f, void* __restrict p){
	printf("Integrand\n");
	T *integrand=(T*)p;
	double  r=intv[0];
	printf("%.3e\n",r);
	getchar();
	double  res=0;
	res+=integrand[0](r);
	*f=res;
	printf("%.3e\n",res);
	return 0;
}
*/

/////////////////////////////////////////////////////
//
////////////////////////////////////////////////////

class testfunc0{
	
	public:
		double operator()(double r, double* p){
			double qs2=pow(2.0e-4/p[1],0.2); 
			
			return(std::cyl_bessel_j(0,r *p[0] )* (1-exp(-r*r *qs2/4))/r);
		}

};
////////////////////////////////////////////////////
class testfunc{
	Collinear_Gluon xg;
	inline double alpha(double mu2 ){
	
			static double b0= ((double)(33 -2*NF))/(12*PI);
			return( 1/(b0* log(mu2/LQCD2)));//LQCD2 lambda_QCD ^2
		}
	
	public:
		double operator()(double r, double* p){
			double q2=pow(r,-2)+1;
			//printf("xg(%.3e %.3e)=%.3e\n", p[1],q2,xg(p[1],q2,1,0.2));
			double qs2=4*PI*PI*alpha(q2)*xg(p[1],q2,1.0,0.2)/(3*85.0); 
			
			return(std::cyl_bessel_j(0,r *p[0] )* (1-exp(-r*r *qs2/4))/r);
		}

};
////////////////////////////////////////////////////////
class testfunc_cheb{
	Chebyshev1D_Collinear_Gluon *xg;
	
	inline double alpha(double mu2 ){
	
			static double b0= ((double)(33 -2*NF))/(12*PI);
			return( 1/(b0* log(mu2/LQCD2)));//LQCD2 lambda_QCD ^2
		}
	
	public:
		testfunc_cheb(Chebyshev1D_Collinear_Gluon& g){
			xg=&g;
		}
		double operator()(double r, double* p){
			double q2=pow(r,-2)+1;
			//printf("xg(%.3e %.3e)=%.3e\n", p[1],q2,xg(p[1],q2,1,0.2));
			double qs2=4*PI*PI*alpha(q2)*(*xg)(p[1],q2)/(3*85.0); 
			//printf("qs2=%.3e r=%.3e x=%.3e, k=%.3e\n", qs2,r,p[1],p[0]);
			double e=r*r *qs2/4;
			
			double val=(e>1.0e-3)?((1-exp(-e))):(e-e*e/2+pow(e,3)/6);
			val/=r;
			
			val*=(r*p[0]>10)?(sqrt(2.0/(PI*r*p[0]))*cos(r*p[0]-PI/4)):(std::cyl_bessel_j(0,r *p[0] ));
			//val*=std::cyl_bessel_j(0,r *p[0] );
			return(val);
		}

};
////////////////////////////////////////////////////////////
// Main
//////////////////////////////////////////////////////////
int main(int argc, char** argv ){
	double par[10]={0};
	CCIntegral cc=CCprepare(64,"dipole",1,5);
	
	par[0]=atof(argv[1]);
	par[1]=atof(argv[2]);
	double imin,imax;
	Kahn acc=Kahn_init(3);
	double val=0,sum=0;
	const double scale=25*(2*PI)/par[0];
	
	double total=0;
	Chebyshev1D_Collinear_Gluon xg;
	
	testfunc_cheb func1(xg);
	testfunc func2;
	testfunc0 func3;

/*
	Kahn_clear(acc);
	sum=0;
	imin=1.0e-8;
	imax=PI/(4*par[0]);
		printf("\n***************************\nGBW\n****************************\n");
	for (int i=0;i<25;++i){
		imax+=scale;
		val=dclenshaw<testfunc0,double*>(cc, func3, par,imin,imax,1.0e-10,1.0e-10);
		sum+=val;
		acc+=val;
		total=Kahn_total(acc);
		printf("[%.3e %.3e ] kt2= %.3e x= %.3e value= %.3e \tsum=%.3e,%.3e %.3e\n" , imin,imax,pow(par[0],2),par[1], val,sum,total,sum-total);
		imin=imax;
	}

	
	sum=0;
	imin=1.0e-8;
	imax=PI/(4*par[0]);
	Levin lev(1500);
		printf("\n***************************\nLevin\n****************************\n");
	for (int i=0;i<50;++i){
		imax+=scale;
		val=dclenshaw<testfunc0,double*>(cc, func3, par,imin,imax,1.0e-10,1.0e-15);
		sum+=val;
		lev.add_term(val);
		
		
		if(i>40){
			printf("[%.3e %.3e ] kt2= %.3e x= %.3e value= %.3e \tsum=%.3e\n" , imin,imax,pow(par[0],2),par[1], val,sum);
			printf("levin= %.3e\n", lev.accel(i-1-10,10));
		}
		imin=imax;
	}
	*/
	/*
	Kahn_clear(acc);
	sum=0;
	imin=1.0e-8;
	imax=10*PI/(4*par[0]);
	printf("\n***************************\nexact\n****************************\n");
	for (int i=0;i<10;++i){
		imax+=scale;
		val=dclenshaw<testfunc,double*>(cc, func2, par,imin,imax,1.0e-10,1.0e-10);
		sum+=val;
		acc+=val;
		total=Kahn_total(acc);
		printf("[%.3e %.3e ] kt2= %.3e x= %.3e value= %.3e \tsum=%.3e,%.3e %.3e\n" , imin,imax,pow(par[0],2),par[1], val,sum,total,sum-total);
		imin=imax;
	}
	*/
	Levin lev(1500);
	//lev.reset();
	Kahn_clear(acc);
	sum=0;
	imin=1.0e-8;
	imax=PI/(4*par[0]);
	xg.init(1.0e-8,1,0.9,pow(imin,-2)*2,1.0,0.2 );
	xg.set_x(par[1]);
	printf("\n***************************\nchebyshev-levin\n****************************\n");
	//FILE* file=fopen("plot.txt","w");
	FILE* file2=fopen("plot2.txt","w");
	for (int i=0;i<50;++i){
		imax+=scale;
		val=dclenshaw<testfunc_cheb,double*>(cc, func1, par,imin,imax,1.0e-15,1.0e-15);
		//val=dclenshaw<testfunc_cheb,double*>(cc, func1, par,(i==0)?(imin):(imax+i*scale),imax+(i+1)*scale,1.0e-15,1.0e-15);
		//printf("%.3e, %.3e\n",func1((imax+imin)/2,par),func1((imax+imin)/2+(imax-imin)/4,par) );
		sum+=val;
		acc+=val;
		total=Kahn_total(acc);
		
		lev.add_term(val);
		if(5*(i/5)==i&&i>15){
			printf("%.3e, %.3e\n",func1((imax+imin)/2,par),func1(imin+(imax-imin)/4,par) );
			printf("[%.3e %.3e ] kt2= %.3e x= %.3e value= %.3e \tsum=%.3e,%.3e %.3e\n" , imin,imax,pow(par[0],2),par[1], val,sum,total,sum-total);
			printf("levin= %.3e\n", lev.accel(i-1-5,5));
		}
		fprintf(file2,"%d %.3e %.3e\n",i,val,sum);
		/*if(i>=10){
			for(int j=0;j<1000;++j){
				fprintf(file,"%.3e %.3e\n",imin+j*(imax-imin)/999,func1(imin+j*(imax-imin)/999,par) );
			}
			
		}*/
		imin=imax;
	}
	//fclose(file);
	fclose(file2);
		
	return 0;
}
	
