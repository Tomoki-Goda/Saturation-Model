
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
	typedef Dipole_Gluon Gluon;
	typedef Sigma SIGMA;

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

	char filenames[500];
	if(opt.path==""){
		opt.path=getenv("DIR");
	}	
	sprintf(filenames,"%s/%s",opt.path.c_str(),opt.input_file_name.c_str());

	FILE* infile=fopen(filenames,"r");
	read_parameters(infile,param);
	fclose(infile);
	double sigpar[10],sudpar[10];


	parameter(param,sigpar,sudpar);

	SIGMA sigma;
	sigma.init(sigpar);
	
	printf("x= %.3e Q2=%.3e\n",opt.x,opt.Q2);	
	sigma.set_kinem(opt.x);
	
	sprintf(filenames,"%s/%s-%.0lf-%.0lf.txt",opt.path.c_str(),"dipole",-log10(opt.x),opt.Q2);
	FILE* outfile=fopen(filenames,"w");
	if(outfile==NULL){
		printf("can't open");
	}else{
		printf("OK: %s\n",filenames);
	}
	double r;
	for(int i =0;i<=100;++i){
		//r=R_MIN*pow(R_MAX/R_MIN, ((double)i)/(100-1));
		r=1.0e-2*pow(30/1.0e-2, ((double)i)/(100-1));
		fprintf(outfile,"%.5e\t%.5e\n",r,sigma(r)/sigpar[0]);

	}
	fclose(outfile);
	
	//printf("\033[1A\033[2K\r");
	sprintf(filenames,"%s/%s-%.0lf-%.0lf.txt",opt.path.c_str(),"gluon", -log10(opt.x),opt.Q2);
	
	outfile=fopen(filenames,"w");
	if(outfile==NULL){
		printf("can't open");
	}else{
		printf("OK: %s\n",filenames);
	}
	Gluon gluon;
	gluon.init(N_APPROX+250,sigpar);
	gluon.set_x(opt.x);
	double k2;
	for(int i =0;i<=100;++i){
		k2=1.0e-2*pow(1.0e+4, ((double)i)/(100-1));
		fprintf(outfile,"%.5e\t%.5e\n",k2,gluon(k2,opt.Q2)/sigpar[0]);

	}
	fclose(outfile);

	return 0;
	
}
