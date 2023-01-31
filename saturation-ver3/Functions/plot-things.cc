#include<iostream>
#include<fstream>
#define TEST 0

#include<cuba.h>
#include<cmath>
#include <vector>
#include <string>
#include"control.h"
#include"./control-default.h"
#include"./constants.h"
#include"Parameters.hh"
#include"./r-formula.hh"
#ifndef GLUON_APPROX
	#define GLUON_APPROX 1
#endif
//#if GLUON_APPROX==1
//	#include"./interpolation.hh"
//	typedef Approx_aF Gluon ;
//	typedef Laplacian_Sigma SIGMA;
//#else
	#include"./gluon-gbw.hh"
	typedef Gluon_GBW Gluon ;
	typedef Sigma SIGMA;

//#endif
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
//#if GLUON_APPROX==0
	sigma.init(sigpar);
//#elif GLUON_APPROX==1
//	sigma.init(N_CHEB_R,sigpar,'s');
//#endif
//
	
	sigma.set_kinem(opt.x);
	
	sprintf(filenames,"%s/%s-%.0lf-%.0lf",opt.path.c_str(),"dipole",opt.x,opt.Q2);
	FILE* outfile=fopen(filenames,"w");
	double r;
	for(int i =0;i<=100;++i){
		r=R_MIN*pow(R_MAX/R_MIN, ((double)i)/(100-1));
		fprintf(outfile,"%.5e\t%.5e\n",r,sigma(r));

	}
	fclose(outfile);
	//printf("\033[1A\033[2K\r");
	sprintf(filenames,"%s/%s-%.0lf-%.0lf",opt.path.c_str(),"gluon",opt.x,opt.Q2);
	outfile=fopen(filenames,"w");
	if(outfile==NULL){
		printf("can't open");
	}else{
		printf("OK: %s\n",opt.output_file_name.c_str());
	}
	Gluon gluon;
//#if GLUON_APPROX==1
//			gluon.init(N_APPROX+50,N_APPROX+50,sigpar);
//			gluon.set_max(5.0e+4);
//#else
			gluon.init(sigpar);
//#endif//GLUON_APPROX==1			

	double k2;
	for(int i =0;i<=100;++i){
		k2=1.0e-7*pow(1.0e+11, ((double)i)/(100-1));
		fprintf(outfile,"%.5e\t%.5e\n",k2,gluon(opt.x,k2,opt.Q2));

	}
	fclose(outfile);

	return 0;
	
}
