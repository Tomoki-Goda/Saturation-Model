//class Approx_aF;
//typedef struct {int i, j; Approx_aF* ptr; } parallel_arg;

/////////////////////////////////////////////////////
// Approx_aF<gluon> af //gluon is Dipole_Gluon or Gluon_GBW
// 
#include<cmath>
#include<iostream>
#include<vector>
#include<string>
#include<complex>
#include<chrono>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_spline.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include<gsl/gsl_dht.h>
#include<gsl/gsl_deriv.h>
#include"control.h"
#include"control-default.h"
#include"constants.h"
#include"clenshaw.hh"
#include"dipole-gluon.hh"

#if GLUON_APPROX==1
typedef Dipole_Gluon GLUON;
#elif MODEL==0
typedef Gluon_GBW GLUON;
#endif


class Approx_aF{
	private:
		GLUON *aF;

		double max_prev=0;
		
		int kt2_npts,x_npts;
		gsl_interp_accel *x_accel_ptr, *kt2_accel_ptr;
		gsl_spline2d *  spline_ptr;
		double *kt2_array=NULL,*x_array=NULL,*aF_array=NULL;
		double mu02=0;
		double sigma_0=0;
		double kt2min=KT2_MIN/2 ,kt2max=-1;
		double xmin=X_MIN;
		const double mu2[1]={0};
		
		int alloc_flag=0;
		
		void free_approx();
		void alloc(int x_npts,int kt2_npts);
		int approximate(const double kt2max);
	
	public:
		double saturation(double x,double kt2_start);
		int export_grid(FILE*file)const;

		Approx_aF(GLUON& g){
			aF=&g;
		}
		~Approx_aF(){
			free_approx();
		}
#if SUDAKOV>=1
		void set_max(double kt2max,const double& mu2);
#else
		void set_max(double kt2max);
#endif
		void init(const int npts1, const int npts2, const double * const &par);
		double operator()(const double x,const double kt2,const double mu2)const;
};
