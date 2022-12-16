#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<vector>
#include"./control.h"
#include"./control-default.h"
#include"./constants.h"
#include"./Parameters.hh"

double INT_PREC=DGAUSS_PREC;
#include"./fcn.h"


int main(int argc, char** argv){

	KtFCN theFCN("/home/tomoki/Saturation-Model/saturation-ver3/data/hera_tot.dat");
	ROOT::Minuit2::MnUserParameters upar;

#if MU202!=0
#if (MODEL==3 && INDEPENDENT_RMAX==0)
	printf(" mu202 is not free. Cannot be controlled\n");
#elif MU0==0
	printf(" mu202 is not used rmax is. mu202 cannot be controlled\n");
#else
	//int pos;
	for(int i=0;i<N_PAR;i++){
		if(strcmp(par_name[i],"mu202")==0){
			par_start[i]=MU202;
			printf("mu202 set to %.5e\n",par_start[i]);
			break;
		}
	}
#endif
#endif
	int skip=0;
	for(unsigned i=0;(i-skip)<N_PAR;i++){
#if MODEL==3
///////////////
//parameters  may be shared bet. BGK & Sud, see dipole-cross-section.h parameters() for how they are organized.
///////////////
#if INDEPENDENT_C==0
		if(i==5){
			skip++;
			continue;
		}
#endif
#if INDEPENDENT_RMAX==0
		if(i==6){
			skip++;
			continue;
		}
#endif
#endif
		upar.Add(par_name[i], par_start[i],par_error[i]);
		upar.SetLimits(par_name[i],par_min[i],par_max[i]);//use migrad.removeLimits(<name>);
	}
	
	//ROOT::Minuit2::MnUserParameterState stat(upar);
	//ROOT::Minuit2::MnStrategy strategy(0);
	//std::cout<<"Parameters\n ";
	//std::cout<<stat<<std::endl;
	INT_PREC=1.0e-3;
	
	ROOT::Minuit2::MnSimplex simplex(theFCN,upar,0);
	ROOT::Minuit2::FunctionMinimum min=simplex(100,10);
	std::cout<<"Parameters "<<min.UserState()<<std::endl;
	
	ROOT::Minuit2::MnMachinePrecision prec;
	prec.SetPrecision(1.0e-5);
	ROOT::Minuit2::MnMigrad migrad(theFCN, min.UserParameters() ,1);
	//ROOT::Minuit2::MnMigrad migrad(theFCN, upar,1);
	INT_PREC=1.0e-5;
	
	//std::cout<<"min= "<<min<<std::endl;
	for(int i=0;i<(N_PAR-skip);i++ ){
		migrad.RemoveLimits(i);
	}
	//ROOT::Minuit2::FunctionMinimum min=migrad();

	//ROOT::Minuit2::FunctionMinimum 
	min=migrad(1000,1);
	
	//stat(upar);
	std::cout<<"Parameters "<<min.UserState()<<std::endl;
	
	std::cout<<"min= "<<min<<std::endl;
	
	//FILE* file=fopen(argv[1],"w");
	std::fstream file;
	file.open(((std::string)argv[1])+"/result.txt",std::fstream::out);
	file<<std::scientific;
	file<<"Qup"<<"\t"<<Q2_MAX<<std::endl;
	for(int i=0;i<(N_PAR-skip);i++){
		//fprintf(file,"%s\t%.10e\n",
		 file<<par_name[i]<<"\t"<<min.UserState().Value(i)<<"\t"<<min.UserState().Error(i)<<std::endl;
	}
	file<<"chisq"<<"\t"<<min.UserState().Fval()<<"\t"<<min.UserState().Edm()<<std::endl;
	file<<"n_data\t"<<theFCN.MAX_N<<std::endl;
	file<<"chisq/dof\t"<<min.Fval()/(theFCN.MAX_N-(N_PAR-skip) ) <<std::endl;
	file.close();
	
	
	
	return 0;
	
}
