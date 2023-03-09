#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<vector>
#include<ctime>
#include<chrono>
#include"control.h"
#include"control-default.h"
#include"constants.h"
#include"Parameters.hh"

int N_APPROX=250;
double INT_PREC=1.0e-4;
#include"kt-formula.hh"
#include"./fcn.h"
#include"../Utilities/options.h"


int main(int argc, char** argv){
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
	printf("********************Program Started********************.\n");
	printf("MODEL     = %d, N_CHEB_R     = %d,             \t IBP          = %d\n",MODEL, N_CHEB_R,IBP);
	printf("R_FORMULA = %d, MU02         = %d\n",R_FORMULA, MU02);
	printf("MODX      = %d, R_CHANGE_VAR = %d,             \t GLUON_APPROX = %d\n",MODX,R_CHANGE_VAR,GLUON_APPROX);
	printf("ALPHA_RUN = %d, Hankel       = %d\n", ALPHA_RUN,  HANKEL);
	printf("FREEZE_QS2= %d, ADD_END      = %d,             \t THRESHOLD    = %d\n",FREEZE_QS2,ADD_END,THRESHOLD);
	printf("Directory = %s ,\t R  = [%.1e, %.1e]\n",opt.path.c_str(),R_MIN,R_MAX);
	printf("*******************************************************.\n");

	KtFCN theFCN("/home/tomoki/Saturation-Model/saturation-ver3/data/hera_tot.dat",opt.path.c_str());
	N_APPROX=250;
	INT_PREC=1.0e-4;
	theFCN.flag=1;
	theFCN(param);
	return 0;	
}
