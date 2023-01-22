#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include"clenshaw.h"


extern double change_var(double & var,double &  jac,const double min, const double max,const double c);
extern double INT_PREC;

//template <typename functype> 
class Interpolate2d{
	protected:
		int x_npts,y_npts;
		gsl_interp_accel *  x_accel_ptr,*  y_accel_ptr;
		gsl_spline2d *  spline_ptr;
		double *x_array,*y_array,*f_array;
		double func(double x ,double y){
			return 0;
		};
		
		void init(int npts1,int npts2 ){
			x_npts=npts1;
			y_npts=npts2;
			x_array=(double*)malloc(x_npts*sizeof(double));
			y_array=(double*)malloc(y_npts*sizeof(double));
			f_array=(double*)malloc(x_npts*y_npts*sizeof(double));
			//func=funcptr;
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
		Interpolate2d(int npts1,int npts2 ){
			init(npts1,npts2);
		}
		
		~Interpolate2d(){
			free_approx();
		}
		
		int approximate(){
			generate_array();
			
			for (int j = 0; j < y_npts; j++){
				for (int i = 0; i < x_npts; i++){
					f_array[j*x_npts+i] = func(x_array[i],y_array[j]);
				}
			}
			x_accel_ptr = gsl_interp_accel_alloc ();
			y_accel_ptr = gsl_interp_accel_alloc ();
			spline_ptr = gsl_spline2d_alloc (gsl_interp2d_bicubic, x_npts,y_npts); // cubic spline
			gsl_spline2d_init (spline_ptr, x_array, y_array, f_array, x_npts, y_npts);
			return(0);
		}
		double approx_f(const double x,const double y)const{	
			double val = gsl_spline2d_eval(spline_ptr, x, y, x_accel_ptr,y_accel_ptr);
			return(val);
		}	
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//               Gluon
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Laplacian_Sigma:public Interpolate2d{
		const double * sigpar;

		double func(const double r,const double x)const{//,const double Q2,const double*sigpar)const {
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
		double approx_f(const double x,const double y)const{	
			double val = gsl_spline2d_eval_deriv_xx(spline_ptr, x, y, x_accel_ptr,y_accel_ptr);
			val += gsl_spline2d_eval_deriv_x(spline_ptr, x, y, x_accel_ptr,y_accel_ptr)/x;
			return(val);
		}
	
	public:
		Laplacian_Sigma(int n_x,int n_y);
		//}
		~Laplacian_Sigma(){//:Interpolate2d(){
		}
		void approx_sigma(){
			approximate();
			
		}
		void set_par(const double * par){
			sigpar=par;
		}
		
		double operator()(const double *R,const std::vector<double>& par)const{
			const double r=*R;
			//const double *par=((const double*)K);
			
			double val=approx_f(r,par[1])*r*2*PI*std::cyl_bessel_j(0,r*sqrt(par[0]));
			//printf("af= %.3e\n",val);
			return( val );
		}
	
};
Laplacian_Sigma::Laplacian_Sigma(int n_x,int n_y): Interpolate2d(n_x,n_y){
}
class Gluon{
	private:
		Laplacian_Sigma dpptr(30,30);
		
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
			double val=3.0/(4*PI*PI)*dclenshaw< Laplacian_Sigma ,std::vector<double>& >(dpptr,par,R_MIN,R_MAX,INT_PREC);
			
			if(std::isnan(val)==1){
				return(0);
			}
#if ALPHA_RUN==1
			val*=alpha(mu2)/0.2;//alpha at mu=1
#endif
			return (val) ;
		}
};

		
















