#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include"clenshaw.h"


extern double change_var(double & var,double &  jac,const double min, const double max,const double c);
extern double INT_PREC;

template <typename functype> 
class Interpolate2d{
	protected:
		int x_npts,y_npts;
		gsl_interp_accel *  x_accel_ptr,*  y_accel_ptr;
		gsl_spline2d *  spline_ptr;
		double *x_array,*y_array,*f_array;
		functype* funcptr;
		
		void init(int npts1,int npts2,const functype &functor ){
			x_npts=npts1;
			y_npts=npts2;
			x_array=(double*)malloc(x_npts*sizeof(double));
			y_array=(double*)malloc(y_npts*sizeof(double));
			f_array=(double*)malloc(x_npts*y_npts*sizeof(double));
			funcptr=&functor;
		}
		void free_approx(){
			gsl_spline2d_free (spline_ptr);
			gsl_interp_accel_free (x_accel_ptr);
			gsl_interp_accel_free (y_accel_ptr);
			free(x_array);
			free(y_array);
			free(f_array);
		}
		
		int generate_array(){
			double x,y;
			for (int j = 0; j < y_npts; j++){
				y=pow(10,-10+10*((double)j)/(y_npts-1));
				y_array[j] = y;
			}	
			for (int i = 0; i < x_npts; i++){
				//x=((double)i)/(x_npts-1);
				x=pow(10,-10+10*((double)i)/(x_npts-1));
				x_array[i] = x;
			}
			return 0;
		}
				
	public:
		Interpolate2d(int npts1,int npts2,const functype &functor ){
			init(npts1,npts2,functor);
		}
		
		~Interpolate2d(){
			free_approx();
		}
		
		int approximate(){
			generate_array();
			
			for (int j = 0; j < y_npts; j++){
				for (int i = 0; i < x_npts; i++){
					f_array[j*x_npts+i] = (*funcptr)(x_array[i],y_array[j]);
				}
			}
			x_accel_ptr = gsl_interp_accel_alloc ();
			y_accel_ptr = gsl_interp_accel_alloc ();
			spline_ptr = gsl_spline2d_alloc (gsl_interp2d_bicubic, x_npts,y_npts); // cubic spline
			gsl_spline2d_init (spline_ptr, x_array, y_array, f_array, x_npts, y_npts);
			return(0);
		}
		double operator()(const double x,const double y)const{	
			double val = gsl_spline2d_eval(spline_ptr, x, y, x_accel_ptr,y_accel_ptr);
			return(val);
		}	
};


template <typename functype> 
class Interpolate_sigma:public Interpolate2d<functype>{
		//using Interpolate2d<functype>::Interpolate2d;
		int generate_array(){
			double x,y;
			for (int j = 0; j < Interpolate2d<functype>::y_npts; j++){
				y=pow(10,-10+10*((double)j)/(Interpolate2d<functype>::y_npts-1));
				Interpolate2d<functype>::y_array[j] = y;
			}	
			for (int i = 0; i < Interpolate2d<functype>::x_npts; i++){
				//x=((double)i)/(x_npts-1);
				x=pow(10,-10+10*((double)i)/(Interpolate2d<functype>::x_npts-1));
				Interpolate2d<functype>::x_array[i] = x;
			}
			return 0;
		}
		public:
		//using Interpolate2d<functype>::Interpolate2d;
		Interpolate_sigma(int npts1,int npts2,const functype &functor ):Interpolate2d<functype>::Interpolate2d(npts1,npts2,functor){
			//Interpolate2d<functype>::init(npts1,npts2,functor);
		}
		
		~Interpolate_sigma(){
			Interpolate2d<functype>::free_approx();
		}
		
		double operator()(const double x,const double y)const{	
			double val = gsl_spline2d_eval_deriv_xx(Interpolate2d<functype>::spline_ptr, x, y, Interpolate2d<functype>::x_accel_ptr,Interpolate2d<functype>::y_accel_ptr);
			val += gsl_spline2d_eval_deriv_x(Interpolate2d<functype>::spline_ptr, x, y, Interpolate2d<functype>::x_accel_ptr,Interpolate2d<functype>::y_accel_ptr)/x;
			return(val);
		}
};

template <typename functype> 
class Interpolate_af:public Interpolate2d<functype>{
		using Interpolate2d<functype>::Interpolate2d;
		int generate_array(){
			double x,y;
			for (int j = 0; j < Interpolate2d<functype>::y_npts; j++){
				y=pow(10,-10+10*((double)j)/(Interpolate2d<functype>::y_npts-1));
				Interpolate2d<functype>::y_array[j] = y;
			}	
			for (int i = 0; i < Interpolate2d<functype>::x_npts; i++){
				x=pow(10,-10+10*((double)i)/(Interpolate2d<functype>::x_npts-1));
				Interpolate2d<functype>::x_array[i] = x;
			}
			return 0;
		}
		public:
		
		
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//               Gluon
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Dipole{
		const double * sigpar;
	public:
		void set_par(const double * par){
			sigpar=par;
		}
		double operator()(const double r,const double x)const{//,const double Q2,const double*sigpar)const {
			double sigma_0=sigpar[0];
			double lambda=sigpar[1];
			double x_0=sigpar[2];
			if(x_0<1.0e-5||x_0>1.0e-3){
				return 0;
			}
			if(lambda<0.05||lambda>0.95){
				return 0;
			}
			
			return( sigma_0*(1-exp( - pow(r * Q0, 2) * pow(x_0/x, lambda)/4)) );	
		}


};
class Approx_Dipole{
	interpolate_sigma<>(30,30,)
	public:
		void approx_sigma(){
			laplacian_sigma.approximate();
		}
		void set_par(const double * par){
			sigpar=par;
		}
		
		double operator()(const double *R,const std::vector<double>& par)const{
			const double r=*R;
			//const double *par=((const double*)K);
			
			double val=lap_sigma(r,par[1])*r*2*PI*std::cyl_bessel_j(0,r*sqrt(par[0]));
			//printf("af= %.3e\n",val);
			return( val );
		}
	
};

class Gluon{
	private:
		const double * sigpar; 
		double mu02=MU02;
		
	public:
		std::string key;
		explicit Gluon(std::string &type,const double (&par)[]){
				key=type;	
				set_par(par);	
		}
		void set_par(const double*par){
			if(key=="gbw"){
#if MU02==0
				mu02=par[3];
#endif
				sigpar=par;
				dpptr.set_par(par);
			}else{
				std::cout<<"unknown model: "<<key<<std::endl;
			}

		}
		~Gluon(){
			
		}
		int set_kinem(const double &b){
			//dpptr.approx_sigma(b);
			dpptr.approx_sigma();
			Q2=b;
			return 0;
		}
	public:
		inline double alpha(const double mu2)const{
			//const double mu02=MU02;
			//return 4.0/(9.0 *std::log( (mu2>mu0)?(mu2/LQCD2):(mu0/LQCD2) ));
			
			return 4.0/(9.0 *std::log( std::max(mu2,mu02)/LQCD2) );
		}
		
	public:
		double operator()(const double x,const double k2,const double mu2)const{
			if(mu02<2*LQCD2){
				return(0);
			}
			std::vector<double> par={k2,x};
			double val=3.0/(4*PI*PI)*dclenshaw< Dipole ,std::vector<double>& >(dpptr,par,R_MIN,R_MAX,INT_PREC);
			
			if(std::isnan(val)==1){
				return(0);
			}
#if ALPHA_RUN==1
			val*=alpha(mu2)/0.2;//alpha at mu=1
#endif
			return (val) ;
		}
};

		
















