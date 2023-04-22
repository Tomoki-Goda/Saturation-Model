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
#include"gluon-integrand.hh"
#include"dipole-gluon.hh"
#include"interpolation-gluon.hh"
#include"options.h"
#include"miscellaneous.hh"



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

	
		SIGMA sigma;
		sigma.init(sigpar);
		INTEG dsigma(sigma);
		dsigma.init(sigpar,'l');
 		

		GLUON dipole_gluon(dsigma);
		dipole_gluon.init(sigpar);
		Approx_aF gluon(dipole_gluon);
		gluon.init(150,150,sigpar);
		const double kt2max=1.0e+6;
	 	printf("Initialized\n");
		gluon.set_max(kt2max);
		printf("Export\n");	
        if(opt.output_file_name!="out.txt"){
            sprintf(outfilenames,"%s",opt.output_file_name.c_str());
        }else{
		    sprintf(outfilenames,"%s/dipole-gluon.txt",opt.path.c_str());
		}
		printf("%s\n",outfilenames);
  		FILE *outfile=fopen(outfilenames,"w");
		gluon.export_grid(outfile);
		printf("Gluon end\n");
		fclose(outfile);
		
		double val,x,r,kt2;
		sprintf(outfilenames,"%s/saturation.txt",opt.path.c_str());
		outfile=fopen(outfilenames,"w");
		for(int i=0;i<10;i++){
			x=1.0e-8*pow(0.01/1.0e-8,((double)i)/9);
			val=gluon.saturation(x,0.1);
			fprintf(outfile,"%.5e\t%.5e\n",x,val);
		}
		fclose(outfile);
		printf("Saturation end\n");
		
		int plot_pts=100;
/////////GLUON///////////////////
		sprintf(outfilenames,"%s/gluon-2.txt",opt.path.c_str());
		outfile=fopen(outfilenames,"w");
		x=1.0e-2;
		for(int i=0;i<plot_pts;i++){
			kt2=1.0e-3*pow(1.0e+6,((double)i)/(plot_pts-1));
			val=gluon(x,kt2,0);
			fprintf(outfile,"%.5e\t%.5e\n",kt2,val/sigpar[0]);
		}
		fclose(outfile);
		sprintf(outfilenames,"%s/gluon-6.txt",opt.path.c_str());
		outfile=fopen(outfilenames,"w");
		x=1.0e-6;
		for(int i=0;i<plot_pts;i++){
			kt2=1.0e-3*pow(1.0e+6,((double)i)/(plot_pts-1));
			val=gluon(x,kt2,0);
			fprintf(outfile,"%.5e\t%.5e\n",kt2,val/sigpar[0]);
		}
		fclose(outfile);
		printf("Gluon Plot end\n");
//////////SIGMA/////////////////////////		
		sprintf(outfilenames,"%s/sigma-2.txt",opt.path.c_str());
		outfile=fopen(outfilenames,"w");
		x=1.0e-2;
		sigma.set_x(x);
		for(int i=0;i<plot_pts;i++){
			r=1.0e-6*pow(5.0e+6,((double)i)/(plot_pts-1));
			val=sigma(x,r);
			fprintf(outfile,"%.5e\t%.5e\n",x,val/sigpar[0]);
		}
		fclose(outfile);
		sprintf(outfilenames,"%s/sigma-6.txt",opt.path.c_str());
		outfile=fopen(outfilenames,"w");
		x=1.0e-6;
		sigma.set_x(x);
		for(int i=0;i<plot_pts;i++){
			r=1.0e-6*pow(5.0e+6,((double)i)/(plot_pts-1));
			val=sigma(x,r);
			fprintf(outfile,"%.5e\t%.5e\n",x,val/sigpar[0]);
		}
		fclose(outfile);
		printf("SIGMA Plot end\n");
	return 0;
	
}
