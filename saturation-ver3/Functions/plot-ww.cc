#include<iostream>
#include<fstream>
#include<cuba.h>
#include<cmath>
#include <vector>
#include <string>
#include"control.h"
#include"control-default.h"
#include"constants.h"
#include"Parameters.hh"

#include<chrono>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_sf.h>
#include<pthread.h>
#include<cuba.h>

extern double  INT_PREC;
#include"clenshaw.hh"
#include"gauss.hh"

#include"gluons.hh"
#include"r-formula.hh"
#include"interpolation-dipole.hh"
#include"dipole-gluon.hh"
#include"interpolation-gluon.hh"
#include "types.hh"
#include"options.h"
//#include"../Utilities/plot.c"


int N_APPROX=N_CHEB_R;

double INT_PREC=1.0e-5;

int main(int argc , char** argv){
	options opt=read_options(argc, argv);
	std::vector<double> param(10,0);

	char infilenames[500];
	char filenames[500];
	if(opt.path==""){
		opt.path=getenv("DIR");
	}
	if(opt.input_file_name!="result.txt"){	
		sprintf(infilenames,"%s",opt.input_file_name.c_str());
	}else{
		sprintf(infilenames,"%s/%s",opt.path.c_str(),opt.input_file_name.c_str());
	}
	
	printf("%s\n",infilenames);
	FILE* infile=fopen(infilenames,"r");
	read_parameters(infile,param);
	fclose(infile);
	double sigpar[10]={0},sudpar[10]={0};
	parameter(param,sigpar,sudpar);
////////////////////////////////////////////////////////////////////////////	
		#if SIGMA_APPROX==-1||SIGMA_APPROX==0//if sigma is not in 1d grid
			SIGMA sigma;
			sigma.init(sigpar);
			
#if SUDAKOV>=1
			Sudakov sudakov;
			sudakov.init(sudpar);
			DSIGMA dsigma(sigma,sudakov);
#else
			DSIGMA dsigma(sigma);
#endif
			dsigma.init(sigpar,'w');
 		#else //-2 and 1  
 			SIGMA sigma;
			sigma.init(sigpar);
#if SUDAKOV>=1
			Sudakov sudakov;
			sudakov.init(sudpar);
			DSIGMA dsigma(sigma,sudakov);
#else
			DSIGMA dsigma(sigma);
#endif
			dsigma.init(N_APPROX+250,sigpar,'w');
		#endif
		
		#if GLUON_APPROX==1
		GLUON dipole_gluon(dsigma);
		dipole_gluon.init(sigpar);
		#else
		GLUON dipole_gluon;
		dipole_gluon.init(sigpar);
		#endif
		Approx_aF<GLUON> gluon(dipole_gluon);
		gluon.init(120,120,sigpar);
		const double kt2max=1.0e+6;
	 	printf("Initialized\n");
	 	//printf("Export\n");	
        if(opt.output_file_name!="out.txt"){
            sprintf(filenames,"%s",opt.output_file_name.c_str());
        }else{
		    sprintf(filenames,"%s/ww.txt",opt.path.c_str());
		}
		printf("%s\n",infilenames);
  		//FILE *outfile=fopen(filenames,"w");
#if SUDAKOV>=1
		FILE *outfile=fopen(filenames,"w");
		fclose(outfile);
		double mu2;
		for(int i=0;i<70;++i){
			mu2=1.0e-1*pow(1.0e+5/1.0e-1,((double)i)/69);
			gluon.set_max(kt2max,mu2);	
			
			outfile=fopen(filenames,"a");
			gluon.export_grid(outfile);
			fclose(outfile);
		}
		printf("Gluon end\n");
#else
		gluon.set_max(kt2max);
  		FILE *outfile=fopen(filenames,"w");
		gluon.export_grid(outfile);
		printf("Gluon end\n");
	 	fclose(outfile);
#endif		

	return 0;
}
