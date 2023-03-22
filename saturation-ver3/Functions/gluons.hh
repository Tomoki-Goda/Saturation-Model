  
//  Claculating of xgpdf(x,Q^2) by DGLAP evolution of initial condition 
//  xgpdf(x,Q^2_0) 
// Originally written by S. Sapeta 
// Modified by T. Goda
// 
// COLGLU Defined at the end
 
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
#include"chebyshev.hh"
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


//extern int N_APPROX;
class Chebyshev_Collinear_Gluon{
	private:
		Collinear_Gluon xg;
	  	cheby cheb[2];
		double A_g=0,l_g=0;
		//int xlen=150,q2len=150;
		//double *xarr=NULL,*q2arr=NULL,*xgarr=NULL;
		double xmin,xmax,q2min,q2max;
		
		double x;
		
	public:
		Chebyshev_Collinear_Gluon(){
			//printf("PrepareChebyshev\n");
			const unsigned deg[]={25,25};
			cheb[0]=PrepareChebyshev(deg,2);
			//printf("PrepareChebyshev End\n");
			
		}
		~Chebyshev_Collinear_Gluon(){
			FreeChebyshev(cheb[0]);
			//FreeChebyshev(cheb[1]);
		}
		double operator()(double* arg,const Collinear_Gluon& xg)const{
			double x=change_var_revert_log(xmin,xmax, arg[0]);//xmin*pow(xmax/xmin,arg[0]);
			double q2=change_var_revert_log(q2min,q2max, arg[1]);//q2min*pow(q2max/q2min,arg[1]);
			return( xg(x,q2,A_g,l_g) );	
		}

		int init(double xmin,double xmax,double q2min, double q2max, double A_g,double l_g ){
			this->xmin=xmin;
			this->xmax=xmax;
			this->q2min=q2min;
			this->q2max=q2max;
			this->A_g=A_g;
			this->l_g=l_g;
			//printf("Chebyshev\n");
			cheb_coeff<Chebyshev_Collinear_Gluon,const Collinear_Gluon>(cheb[0],*this,xg);
			//printf("Chebyshev done\n");
			return 0;
		}
		int set_x(double x){
			this->x=x;
			//cheb_coeff<Chebyshev_Collinear_Gluon,const Collinear_Gluon>(cheb[1],*this,xg);
			//printf("set x\n");
			cheb[1]=chebyshev_reduce(cheb[0], change_var_compactify_log(xmin,xmax,x ), 0 );
			//printf("x set\n");
			return 0;
		}
		double operator()(const double x,const double Q2)const{
			if(x>xmax){
				printf("x too large: %.3e < %.3e < %.3e\n",xmin,x,xmax );
			}
			if(x<xmin){
				printf("x too small: %.3e < %.3e < %.3e\n",xmin,x,xmax );
			}
			if(Q2>q2max){
				printf("Q2 too large: %.3e < %.3e < %.3e\n",q2min,Q2,q2max );
			}
			if(Q2<q2min){
				printf("Q2 too small: %.3e < %.3e < %.3e\n",q2min,Q2,q2max );
			}
#if GLUON_APPROX!=0	///////////////////////////////////////////////////////////////////		
// This block should be removed if x varies often.
// Only useful if x is fixed and Q changes rapidly.
///////////////////////////////////////////////////////////////////////////////////////////
			static double x0=0;
			static cheby chebq2[1];
			double res1;
			if(x0!=x){
				chebq2[0]=chebyshev_reduce(cheb[0], change_var_compactify_log(xmin,xmax,x ), 0 );
				/*double arg[]={
					change_var_compactify_log(xmin,xmax,x ),
					change_var_compactify_log(q2min,q2max,Q2 )
				};
				res0=chebyshev(cheb[0],arg);*/
				x0=x;
			}
			double arg[]={
				change_var_compactify_log(q2min,q2max,Q2 )
			};
			res1=chebyshev(chebq2[0],arg);
			/*if(res1!=res0){
				printf("Chebyshev failed: %.3e %.3e diff= %.3e\n ",res0,res1,res1-res0);
				getchar();
			}*/
			return(res1);
#else //////////////////////////////////////////////////////////////////////////////////////
			double arg[]={
				change_var_compactify_log(xmin,xmax,x ),
				change_var_compactify_log(q2min,q2max,Q2 )
			};
			return(chebyshev(cheb[0],arg));
#endif
			
		}
		
		/*double operator()(const double Q2)const{
			if(Q2>q2max){
				printf("Q2 too large: %.3e < %.3e < %.3e\n",q2min,Q2,q2max );
			}
			if(Q2<q2min){
				printf("Q2 too small: %.3e < %.3e < %.3e\n",q2min,Q2,q2max );
			}
			//double x1=change_var_compactify_log(xmin,xmax,x );
			//double x2=change_var_compactify_log(q2min,q2max,Q2 );
			double arg[1]={
				change_var_compactify_log(q2min,q2max,Q2 )
			};
			
			return(chebyshev(cheb[1],arg));
		}
		*/

};

extern int N_APPROX;
class Interpolate_Collinear_Gluon{
	private:
		Collinear_Gluon xg;
		gsl_interp_accel *x_accel_ptr=NULL, *q2_accel_ptr=NULL;
		gsl_spline2d *  spline_ptr=NULL;
		double A_g=0,l_g=0;
		int xlen=150,q2len=150;
		double *xarr=NULL,*q2arr=NULL,*xgarr=NULL;
		double xmin,xmax,q2min,q2max;
	public:
		//explicit Interpolate_Collinear_Gluon(int xlen,int q2len){
		explicit Interpolate_Collinear_Gluon(){
			//this->q2len=q2len;
			//this->xlen=xlen;

			xarr=(double*)malloc(xlen*sizeof(double));
			q2arr=(double*)malloc(q2len*sizeof(double));
			xgarr=(double*)malloc(xlen*q2len*sizeof(double));
			x_accel_ptr = gsl_interp_accel_alloc ();
			q2_accel_ptr = gsl_interp_accel_alloc ();
			spline_ptr = gsl_spline2d_alloc(gsl_interp2d_bicubic,q2len,xlen);
		}
		~Interpolate_Collinear_Gluon(){
			gsl_spline2d_free (spline_ptr);
			gsl_interp_accel_free (x_accel_ptr);
			gsl_interp_accel_free (q2_accel_ptr);
			free(xarr);
			free(q2arr);
			free(xgarr);
		}
		int init(double xmin,double xmax,double q2min, double q2max, double A_g,double l_g ){
			//this->A_g=A_g;
			//this->l_g=l_g;
#pragma omp parallel
{
//			double x,q2;
#pragma omp for
			for( int i=0;i<xlen;++i){
				double x=xmin*pow(xmax/xmin,((double)i)/(xlen-1));
				xarr[i]=x;
				for(int j=0;j<q2len;++j){
					double q2=q2min*pow(q2max/q2min,((double)j)/(q2len-1));
					if(i==0){
						q2arr[j]=q2;
					}
					xgarr[i*q2len+j]=xg(x,q2,A_g,l_g);
				}
			}
}
			this->xmin=xmin;
			this->xmax=xmax;
			this->q2min=q2min;
			this->q2max=q2max;

			gsl_spline2d_init (spline_ptr,q2arr, xarr, xgarr, q2len, xlen);
			//printf("gluon approxed \n");
			return 0;
		}
		double operator()(const double x,const double Q2)const{
			if(x>xmax){
				printf("x too large: %.3e < %.3e < %.3e\n",xmin,x,xmax );
			}
			if(x<xmin){
				printf("x too small: %.3e < %.3e < %.3e\n",xmin,x,xmax );
			}
			if(Q2>q2max){
				printf("Q2 too large: %.3e < %.3e < %.3e\n",q2min,Q2,q2max );
			}
			if(Q2<q2min){
				printf("Q2 too small: %.3e < %.3e < %.3e\n",q2min,Q2,q2max );
			}
			return(gsl_spline2d_eval(spline_ptr,Q2, x,q2_accel_ptr, x_accel_ptr));
		}

};
//////////////////////////////////////////
//////////////////////////////////////////

#endif
