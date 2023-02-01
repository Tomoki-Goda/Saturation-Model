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
	#include"./interpolation.hh"
	typedef Approx_aF Gluon ;
	typedef Laplacian_Sigma SIGMA;
//#else
//	#include"./gluon-gbw.hh"
// 	typedef Gluon_GBW Gluon ;
//	typedef Sigma SIGMA;

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
double INT_PREC=1.0e-5;
int main(int argc , char** argv){
	options opt=read_options(argc, argv);
	std::vector<double> param(10,0);

	char filenames[500];
	if(opt.path==""){
		opt.path=getenv("DIR");
	}	
	sprintf(filenames,"%s/%s",opt.path.c_str(),opt.input_file_name.c_str());
	printf("%s\n",filenames);
	FILE* infile=fopen(filenames,"r");
	read_parameters(infile,param);
	fclose(infile);
	double sigpar[10],sudpar[10];


	parameter(param,sigpar,sudpar);
	//printf("Start, press any key\n");
	//getchar();
	SIGMA sigma;
	sigma.init(N_CHEB_R+25,sigpar,'s');
	printf("kinem set\n");
	
	sprintf(filenames,"%s/%s",opt.path.c_str(),"dipole-grid.txt");
	FILE* outfile=fopen(filenames,"w");
	double x;
	for(int i=0;i<N_APPROX+50;i++){
		x=pow(10,-15+15*((double)i)/(N_APPROX+50-1));
		sigma.approximate(x);
		sigma.export_grid(outfile);
	}
	fclose(outfile);

	
	sprintf(filenames,"%s/%s",opt.path.c_str(),"gluon-grid.txt");
	outfile=fopen(filenames,"w");
	Gluon gluon;
//#if GLUON_APPROX==1
			gluon.init(N_APPROX+50,N_APPROX+50,sigpar);
			gluon.set_max(5.0e+4);
//#else
//			gluon.init(sigpar);
//#endif//GLUON_APPROX==1			
	gluon.export_grid(outfile);
	fclose(outfile);

	return 0;
	
}
