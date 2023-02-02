/*  
 *  Claculating of xgpdf(x,Q^2) by DGLAP evolution of initial condition 
 *  xgpdf(x,Q^2_0)
 * 
 * Originally written by S. Sapeta 
 * Modified by T. Goda
 */

//#include "main.h"
#include "./complex.hh"
#include "cfortran.h"
#include "./clenshaw.hh"
#include "./gauss.hh"
/* CERNLIB functions*/
extern "C" doublecomplex wgamma_(const doublecomplex*);
extern "C" doublecomplex wpsipg_(const doublecomplex*,int*);
extern "C" double dgammf_(const double*);
class Collinear_Gluon{
	CCIntegral cc;//=CCprepare(64,"gluon",10);
	
	private:
		const double       beta = 6.6;
		const double        n_0 = 0.5;       /* Maximal singluraity of integrand */
		//int flag=0;
		//double A_g, lambda_g;

		doublecomplex gammatilde(const doublecomplex n)const{
			doublecomplex n1,n2,n3,l1,l2,t1,t2,t3,cx,value;
			double m1,m2,rl;
			int k=0;
			n1 = n+1.0;
		    	n2 = n+2.0;
		    	n3 = n+3.0;

		    	l1 = Conjg(n)*Conjg(n1);
		   	m1 = Cabs(n)*Cabs(n)*Cabs(n1)*Cabs(n1);
		   	t1 = (1.0/m1)*l1;

		  	l2 = Conjg(n2)*Conjg(n3);
			m2 = Cabs(n2)*Cabs(n2)*Cabs(n3)*Cabs(n3);
			t2 = (1.0/m2)* l2;

			t3 = wpsipg_(&n2,&k);

			cx = (t1+t2)-t3;
			rl = 11.0/2.0-NF/3.0-6.0*GAMMA_E;

			value = 6.0*cx + rl;
			return value;
		}

	public:
		explicit Collinear_Gluon(){
			cc=CCprepare(128,"gluon",50);
			//printf("gluon\n");
		}
		~Collinear_Gluon(){
			//printf("gluon end\n");
		}
		double operator()(const double y,const double*par )const {
		//doublecomplex xgpdf_integrand(double y, double Y, double t) {
			doublecomplex n0,n1,n2,g1,g2,gt,ex,l;
			double val;
			double m;

			double Yg=par[0], tg=par[1];
			double lambda_g=par[2];

			n0 = Complex(n_0,y);
			n1 = Complex(-lambda_g+n_0,y);
			n2 = Complex(-lambda_g+beta+n_0,y);

			g1 = wgamma_(&n1);
			g2 = wgamma_(&n2);

			gt = tg *gammatilde(n0 );

			ex = Cexp(Complex(0,y* Yg)+gt);

			l = g1*Conjg(g2);
			m = Cabs(g2)*Cabs(g2);

			val = ((1.0/m)*l*ex).r;
			//static int flag=0;
			if(not(std::isfinite(val))){
				//if(flag<=100||(flag/100)*100==flag){
					printf("\033[1A\033[2K\r");
					printf("\033[1A\033[2K\r");
					printf("\033[1A\033[2K\r");
					printf("\033[1A\033[2K\r");
				//}else{
				//	++flag;
				//}
		    		printf("%.3e %.3e  \n",Yg,tg);
		    		printf("%.3e %.3e %.3e %.3e \t %.3e %.3e %.3e \t %.3e   \n",n0.r,n1.r,n2.r, g1.r,g2.r,gt.r,ex.r,l.r);
		    		printf("%.3e %.3e %.3e %.3e \t %.3e %.3e %.3e \t %.3e  %.3e \nval=%.3e\n",n0.i,n1.i,n2.i, g1.i,g2.i,gt.i,ex.i,l.i,m,val);
		    		//getchar();
		    		return 0;
		    	}
		    	
			return val;
		}


		//void set_xg_parameter(double ag,double lg){
		//	A_g=ag;
		//	lambda_g=lg;
		//}

		/*******************************************************************************
		* Gluons pdf  function
		*******************************************************************************/
		double operator()(const double x, const double QQ,const double A_g,const double l_g)const  {
			//return(A_g*pow(x,-l_g));
			static int flag_nan=0;
			double normalization;
			double value;
			const double bprim = 33.0/6.0-NF/3.0;
			double par[] = {
				log(1/x),
				(1/bprim)*log(log(QQ/LQCD2)/log(Q0/LQCD2)),
				l_g
			};
		    	normalization = A_g*exp(n_0* par[0] )*dgammf_(&beta)/PI;
			//value=dclenshaw<const Collinear_Gluon, const double*>(*this,par, a,c,NRel,1.0e-15);
			//value=dgauss<const Collinear_Gluon, const double*>(*this,par, a,c,NRel,1.0e-15); 
			value=dclenshaw< const Collinear_Gluon, const double*>(cc,*this,par, 0,150,1.0e-15,1.0e-17);  
			//this->flag=0;	
			value=normalization*value;
			if(!std::isfinite(value)||value<0){
				if(flag_nan==0){
					std::cout<<std::scientific<<"gluon error:: "<<value<<" for x="<<x<<" Q2= "<<QQ<<std::endl;
					std::cout<<A_g<<"  "<<l_g<<std::endl;
					flag_nan=1;
				}
				return 0;
			}
			//printf("x=%.3e mu2= %.3e Ag= %.3e lg= %.3e\n",x,QQ,A_g,l_g);
			//getchar();
			return value ;
		}

};

