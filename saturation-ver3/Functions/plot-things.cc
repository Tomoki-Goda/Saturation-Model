#include<iostream>
#include<fstream>
#define TEST 0

#include<cuba.h>
#include<cmath>
#include <vector>
#include <string>
#include"./control.h"
#include"./control-default.h"
#include"./constants.h"
#include"Parameters.hh"
#include"./r-formula.hh"
#ifndef GLUON_APPROX
	#define GLUON_APPROX 1
#endif
#if GLUON_APPROX==1
	#include"./interpolation.hh"
	typedef Approx_aF Gluon ;
	typedef Laplacian_Sigma SIGMA;
#else
	#include"./gluon-gbw.hh"
	typedef Gluon_GBW Gluon ;
	typedef Sigma SIGMA;

#endif
//#include"./clenshaw.h"


//#include"./dgauss.h"
//#include"./clenshaw-curtis.hh"
extern PREC INT_PREC;
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

int main(int argc , char** argv){
	options opt=read_options(argc, argv);
	std::vector<double> param(10,0);
	FILE* infile=fopen(opt.input_file_name,"r");
	read_parameters(infile,param);
	fclose(infile);
	double sigpar[10],sudpar[10];


	parameter(param,sigpar,sudpar);
	SIGMA sigma;
#if GLUON_APPROX==0
	sigma.init(sigpar);
#elif GLUON_APPROX==1
	sigma.init(50,sigpar);
#endif

	sigma.set_kinem(opt.x);
	FILE* outfile=fopen(opt.output_file_name,"w");
	double r;
	for(int i =0;i<=100;++i){
		r=R_MIN*pow(R_MAX/R_MIN, ((double)i)/(100-1));
		fprintf(outfile,"%.5e\t%.5e\n",r,sigma(r));

	}
	fclose(outfile);
	return 0;
	
}
