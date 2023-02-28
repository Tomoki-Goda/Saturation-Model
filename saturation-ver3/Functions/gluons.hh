/*  
 *  Claculating of xgpdf(x,Q^2) by DGLAP evolution of initial condition 
 *  xgpdf(x,Q^2_0)
 * 
 * Originally written by S. Sapeta 
 * Modified by T. Goda
 */
#ifndef GLUONS_HH
#define GLUONS_HH
//#include "main.h"

#include<cmath>
//#include "./complex.hh"
//#include "cfortran.h"
#include"./control-default.h"
#include"./constants.h"
#include "./clenshaw.hh"
#include "./gauss.hh"
#include"polygamma.hh"
#include<complex>
#include<pthread.h>
#include<gsl/gsl_sf_gamma.h>
/* CERNLIB functions*/
//extern "C" doublecomplex wgamma_(const doublecomplex*);
//extern "C" doublecomplex wpsipg_(const doublecomplex*,int*);
//extern "C" double dgammf_(const double*);
//pthread_mutex_t mut1, mut2,mut3;


class Collinear_Gluon{
	CCIntegral cc=CCprepare(256,"gluon",1,3);
	private:
		const double       beta = 6.6;
		double dgammafbeta;
		const double        n_0 = 0.5;       /* Maximal singluraity of integrand */

		std::complex<double> gammatilde(const std::complex<double>& n)const{
			std::complex<double> n1,n2,n3,l1,l2,t1,t2,t3,cx,value;
			double m1,m2,rl;
			int k=0;
			n1 = n+1.0;
		    	n2 = n+2.0;
		    	n3 = n+3.0;

		    	l1 = conj(n)*conj(n1);
		   	m1 = abs(n)*abs(n)*abs(n1)*abs(n1);
		   	t1 = (1.0/m1)*l1;

		  	l2 = conj(n2)*conj(n3);
			m2 = abs(n2)*abs(n2)*abs(n3)*abs(n3);
			t2 = (1.0/m2)* l2;

			t3 = digamma(n2);

			cx = (t1+t2)-t3;
			rl = 11.0/2.0-NF/3.0-6.0*GAMMA_E;
			
			value = 6.0*cx + rl;
			return value;
		}

	public:
		//Collinear_Gluon operator=(Collinear_Gluon Collinear_Gluon )
		Collinear_Gluon(const Collinear_Gluon& init){
			dgammafbeta=init.dgammafbeta;
		}
		
		explicit Collinear_Gluon(){
			dgammafbeta=gsl_sf_gamma(beta)/PI;
		}
		~Collinear_Gluon(){
		}
		double operator()(const double y,const std::vector<double> &par)const {
			std::complex<double> n0,n1,n2,g1,g2,gt,ex,l;
			double val;
			double m;

			const double Yg=par[0], tg=par[1];
			const double lambda_g=par[2];
			gsl_sf_result resr,resi;
			std::complex<double>comp(0,1);
			n0 = n_0+y*comp;
			n1 =-lambda_g+n_0+comp*y;
			n2 =-lambda_g+beta+n_0+comp*y;

			gsl_sf_lngamma_complex_e(n1.real(),n1.imag(),&resr,&resi );
			g1=exp(resr.val+comp*resi.val);
			gsl_sf_lngamma_complex_e(n2.real(),n2.imag(),&resr,&resi );
			g2=exp(resr.val+comp*resi.val);
			
			gt = tg *gammatilde(n0 );
			ex = exp(comp*y* Yg+gt);
			l = g1*conj(g2);
			m = abs(g2)*abs(g2);
			val = ((1.0/m)*l*ex).real();
			if(not(std::isfinite(val))){
		    		return 0;
		    	}
			return val;
		}

		/*******************************************************************************
		* Gluons pdf  function
		*******************************************************************************/
		double operator()(const double x, const double QQ,const double A_g,const double l_g)const  {
			static int flag_nan=0;
			double normalization;
			double value;
			const double bprim = 33.0/6.0-NF/3.0;
			const std::vector<double>par{
				log(1/x),
				(1/bprim)*log(log(QQ/LQCD2)/log(Q0/LQCD2)),
				l_g
			};
		    	normalization = A_g*exp(n_0* par[0] )*dgammafbeta;
			value=dclenshaw<const Collinear_Gluon, const std::vector<double> >(cc,*this,par,0,150,1.0e-10,1.0e-15);  
			//this->flag=0;	
			value=normalization*value;
		 	if(!std::isfinite(value)||value<0){
				//if(flag_nan==0){
				//	std::cout<<std::scientific<<"gluon error:: "<<value<<" for x="<<x<<" Q2= "<<QQ<<std::endl;
				//	std::cout<<A_g<<"  "<<l_g<<std::endl;
				//	flag_nan=1;
				//}
				return 0;
			}
			//printf("x=%.3e mu2= %.3e Ag= %.3e lg= %.3e\n",x,QQ,A_g,l_g);
			//getchar();
			//pthread_mutex_unlock(&mut);
			return value ;
		}

};
#endif

