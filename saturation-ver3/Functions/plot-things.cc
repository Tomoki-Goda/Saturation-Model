
#include<iostream>
#include<fstream>
#define TEST 0

#include<cuba.h>
#include<cmath>
#include <vector>
#include <string>
#include"control.h"
#include"control-default.h"
#include"constants.h"
#include"Parameters.hh"
#ifndef GLUON_APPROX
	#define GLUON_APPROX 1
#endif

#include<chrono>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

extern double  INT_PREC;
#include<pthread.h>
#include<cuba.h>

#include"clenshaw.hh"
#include"gauss.hh"

#include"gluons.hh"
#include"r-formula.hh"
#include"interpolation-dipole.hh"
#include"dipole-gluon.hh"
#include"interpolation-gluon.hh"
//#if GLUON_APPROX==1
//	#include"interpolation.hh"
//	typedef Approx_aF Gluon ;
//	typedef Laplacian_Sigma SIGMA;
//#else
//	#include"./gluon-gbw.hh"
//	typedef Gluon_GBW Gluon ;
#if SIGMA_APPROX<0
	#if CHEB_D==2
		typedef Chebyshev_Collinear_Gluon COLGLU;
	#elif CHEB_D==1
		typedef Chebyshev1D_Collinear_Gluon COLGLU;
	#endif
#else 	
		typedef Collinear_Gluon COLGLU;
#endif

#if MODEL==1
	typedef Sigma<COLGLU> SIGMA ;
#else 
	typedef Sigma SIGMA ;
#endif

#if ((SIGMA_APPROX==0)||(SIGMA_APPROX==-1))//negative means xg is approximated.
		typedef Gluon_Integrand<SIGMA> DSIGMA;
#elif ((SIGMA_APPROX>0)||(SIGMA_APPROX<-1))
	typedef Laplacian_Sigma<SIGMA> DSIGMA;
#endif


typedef Dipole_Gluon<DSIGMA> GLUON;

//#endif
//#include"./clenshaw.h"


//#include"./dgauss.h"
//#include"./clenshaw-curtis.hh"
//#include"./r-formula.hh"
//int CUBACORES=4;
//#include"./Photon.hh"

#include"../Utilities/options.h"
//#include"../Utilities/plot.c"


#ifndef ALPHA_RUN
	#define ALPHA_RUN 0 
#endif
#ifndef MODX
	#define MODX 0
#endif
#ifndef PHI
	#define PHI 0
#endif

#ifndef SCATTER
	#define SCATTER 0
#endif

#ifndef MU02
	#define MU02 1
#endif
int N_APPROX=N_CHEB_R;
double INT_PREC=1.0e-4;
int main(int argc , char** argv){
	options opt=read_options(argc, argv);
	std::vector<double> param(10,0);

	char outfilenames[500], infilenames[500];
	if(opt.path==""){
		opt.path=getenv("DIR");
	}	
	sprintf(infilenames,"%s/%s",opt.path.c_str(),opt.input_file_name.c_str());

	FILE* infile=fopen(infilenames,"r");
	read_parameters(infile,param);
	fclose(infile);
	double sigpar[10],sudpar[10];


	parameter(param,sigpar,sudpar);

		#if SIGMA_APPROX==-1||SIGMA_APPROX==0//if sigma is not in 1d grid
			SIGMA sigma;
			sigma.init(sigpar);
			DSIGMA dsigma(sigma);
			dsigma.init(sigpar,'l');
 		#else //-2 and 1  
 			SIGMA sigma;
			sigma.init(sigpar);
			DSIGMA dsigma(sigma);
			dsigma.init(N_APPROX+250,sigpar,'l');
		#endif
		
		GLUON dipole_gluon(dsigma);
		dipole_gluon.init(sigpar);
		Approx_aF<GLUON> gluon(dipole_gluon);
		gluon.init(50,50,sigpar);
		const double kt2max=1.0e+8;
	 	printf("Initialized\n");
		gluon.set_max(kt2max);
		printf("Export\n");	
        if(opt.output_file_name!="out.txt"){
            sprintf(outfilenames,"%s",opt.output_file_name.c_str());
        }else{
		    sprintf(outfilenames,"%s/ww.txt",opt.path.c_str());
		}
		printf("%s\n",outfilenames);
  		FILE *outfile=fopen(outfilenames,"w");
		gluon.export_grid(outfile);
		printf("Gluon end\n");
	 	fclose(outfile);
	return 0;
	
}
