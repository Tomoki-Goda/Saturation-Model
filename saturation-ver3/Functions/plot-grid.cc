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




//#if GLUON_APPROX==1
#if HANKEL==1
	typedef Hankel_aF Gluon ;
#else
	typedef Approx_aF Gluon ;
#endif
	typedef Laplacian_Sigma SIGMA;
//#else
//	#include"./gluon-gbw.hh"
// 	typedef Gluon_GBW Gluon ;
//	typedef Sigma SIGMA;

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
int N_APPROX=250;

double INT_PREC=1.0e-4;

int main(int argc , char** argv){
	options opt=read_options(argc, argv);
	std::vector<double> param(10,0);

	char filenames[500];
	char filenames2[500];
	if(opt.path==""){
		opt.path=getenv("DIR");
	}	
	sprintf(filenames,"%s/%s",opt.path.c_str(),opt.input_file_name.c_str());
	printf("%s\n",filenames);
	FILE* infile=fopen(filenames,"r");
	read_parameters(infile,param);
	fclose(infile);
	double sigpar[10]={0},sudpar[10]={0};
	parameter(param,sigpar,sudpar);
	//printf("Start, press any key\n");
	//getchar();
	SIGMA sigma;
	sigma.init(N_CHEB_R+150,sigpar,'s');
	printf("kinem set\n");
	
	sprintf(filenames,"%s/%s",opt.path.c_str(),"dipole-grid.txt");
	FILE* outfile=fopen(filenames,"w");
	//sprintf(filenames2,"%s/%s",opt.path.c_str(),"laplacian-grid.txt");
	//FILE* outfile2=fopen(filenames2,"w");
	double x;
	for(int i=0;i<N_APPROX;i++){
		x=pow(10,-15+15*((double)i)/(N_APPROX-1));
		sigma.set_kinem(x);
		sigma.export_grid(outfile);
	}
	fclose(outfile);
	//fclose(outfile2);

	
	sprintf(filenames,"%s/%s",opt.path.c_str(),"gluon-grid.txt");
	printf("%s/%s\n",opt.path.c_str(),"gluon-grid.txt");
	
  	outfile=fopen(filenames,"w");
	Gluon gluon;
//#if GLUON_APPROX==1
	printf("Gluon start\n");
	gluon.init(N_APPROX,N_APPROX,N_APPROX+250,sigpar);
	printf("Initialized\n");
	gluon.set_max(5.0e+4);
	printf("Export\n");			
	gluon.export_grid(outfile);
	printf("Gluon end\n");
	fclose(outfile);
	
	double sat;
	sprintf(filenames,"%s/%s",opt.path.c_str(),"saturation.txt");
	printf("%s/%s\n",opt.path.c_str(),"saturation.txt");
  	outfile=fopen(filenames,"w");
	
	for(int i=0;i<10;i++){
		x=pow(10,-8+6*((double)i)/(10-1));
		sat=gluon.saturation(x,0.2);
		fprintf(outfile,"%.5e\t%.5e\n",x,sat );
	}
	fclose(outfile);

	return 0;
	
}
