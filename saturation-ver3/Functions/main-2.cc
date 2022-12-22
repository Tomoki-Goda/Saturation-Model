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

//KtFCN theFCN("/home/tomoki/Saturation-Model/saturation-ver3/data/hera_tot.dat");
int check_min(ROOT::Minuit2::FunctionMinimum *min, int ndata){
	for(int i=0;i<ndata;i++){
		if ((std::isnan(min->UserState().Value(i))+std::isinf(min->UserState().Value(i)) )!=0){
			std::cout<<"Parameter "<<i<<"  "<<par_name[i]<<"  "<<  min->UserState().Value(i) <<std::endl;
			return 1;
		}
	}
	if ((std::isnan(min->UserState().Fval())+std::isinf(min->UserState().Fval()) )!=0){
		printf("Fval = %.3e",min->UserState().Fval() );
		return 1;
	}
	if ((std::isnan(min->UserState().Edm())+std::isinf(min->UserState().Edm()) )!=0){
		printf("Edm = %.3e",min->UserState().Edm() );
		return 1;
	}

	return 0;

}

int save_res(std::string name, const ROOT::Minuit2::FunctionMinimum *min,const KtFCN *theFCN, int ndata){
	std::fstream file;
	file.open(name ,std::fstream::out);
	file<<std::scientific;
	file<<"Qup"<<"\t"<<Q2_MAX<<std::endl;
	for(int i=0;i<ndata;i++){
		//fprintf(file,"%s\t%.10e\n",
		 file<<par_name[i]<<"\t"<<min->UserState().Value(i)<<"\t"<<min->UserState().Error(i)<<std::endl;
	}
	file<<"chisq"<<"\t"<<min->UserState().Fval()<<"\t"<<min->UserState().Edm()<<std::endl;
	file<<"n_data\t"<<theFCN->MAX_N<<std::endl;
	file<<"chisq/dof\t"<<min->Fval()/(theFCN->MAX_N-(ndata) ) <<std::endl;
	file<<"Flag\t"<<min->IsValid()<<std::endl;
	file<<"eps\t"<<INT_PREC<<std::endl;
	
	file.close();
	return 0;
}

int main(int argc, char** argv){
	std::chrono::system_clock walltime;
	std::chrono::time_point start= walltime.now();
	
	std::cout<<std::scientific<<std::endl;
	
	printf("Czesc World!\nProgram Started.\n");
	printf("Directory=%s\n",(char*)argv[1]);
	

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
	
	INT_PREC=1.0e-2;
	printf("*****************************\n");
	printf("*** Simplex: eps=%.1e  ***\n",INT_PREC);
	printf("*****************************\n");
	ROOT::Minuit2::MnSimplex simplex(theFCN,upar,0);
	ROOT::Minuit2::FunctionMinimum min=simplex(200,1);
	ROOT::Minuit2::FunctionMinimum min_prev=min;

	int flag=check_min(&min,N_PAR-skip);
	if( flag==0){
		min_prev=min;
	}else{
		min=min_prev;
	}
	std::cout<<"Parameters "<<min_prev.UserState()<<std::endl;
	std::cout<<"Parameters "<<min.UserState()<<std::endl;
	
	save_res(((std::string)argv[1])+"/result.txt",&min,&theFCN,N_PAR-skip);
	
	INT_PREC=1.0e-3;
	printf("*****************************\n");
	printf("*** Simplex: eps=%.1e  ***\n",INT_PREC);
	printf("*****************************\n");
	min=simplex(100,1);
	
	flag=check_min(&min,N_PAR-skip);
	if( flag==0){
		min_prev=min;
	}else{
		min=min_prev;
	}
	std::cout<<"Parameters "<<min_prev.UserState()<<std::endl;
	std::cout<<"Parameters "<<min.UserState()<<std::endl;
		
	save_res(((std::string)argv[1])+"/result.txt",&min,&theFCN,N_PAR-skip);	

	ROOT::Minuit2::MnMachinePrecision prec;
	prec.SetPrecision(1.0e-8);

	
	ROOT::Minuit2::MnMigrad migrad(theFCN, min.UserParameters() ,0);
	for(int i=0;i<(N_PAR-skip);i++ ){
		migrad.RemoveLimits(i);
	}
	

	INT_PREC=1.0e-3;
	printf("***************************\n");
	printf("*** First: eps=%.1e  ***\n",INT_PREC);
	printf("***************************\n");
	for(int i=0;i<10;i++){
		min=migrad(20,1);
		flag=check_min(&min,N_PAR-skip);
		if( flag==0){
			min_prev=min;
		}else{
			min=min_prev;
		}
		std::cout<<"Parameters "<<min.UserState()<<std::endl;
		std::cout<<std::scientific<<"fcn= "<<min.UserState().Fval()<<", edm= "<<min.UserState().Edm()<<std::endl;
		
		
		if(min.IsValid()){
			printf(" %.3e/%.3e = %.3e\n", min.UserState().Edm(),  (min.UserState().Fval()),min.UserState().Edm()/ (min.UserState().Fval()));
			break;
		}
	}
	
	save_res(((std::string)argv[1])+"/result.txt",&min,&theFCN,N_PAR-skip);	
	INT_PREC=1.0e-4;
	printf("***************************\n");
	printf("*** Second: eps=%.1e  ***\n",INT_PREC);
	printf("***************************\n");
	ROOT::Minuit2::MnMigrad migrad2(theFCN, min.UserParameters() ,1);
	
	for(int i=0;i<20;i++){
		min=migrad2(50,1);
		flag=check_min(&min,N_PAR-skip);
		if( flag==0){
			min_prev=min;
		}else{
			min=min_prev;
		}
		std::cout<<"Parameters "<<min.UserState()<<std::endl;
		std::cout<<std::scientific<<"fcn= "<<min.UserState().Fval()<<", edm= "<<min.UserState().Edm()<<"  "<<min.IsValid()<<std::endl;
		if(min.IsValid()){
			printf("FINE: %.3e/%.3e = %.3e\n", min.UserState().Edm(),  (min.UserState().Fval()),min.UserState().Edm()/ (min.UserState().Fval()));
			break;
		}
	}
	save_res(((std::string)argv[1])+"/result.txt",&min,&theFCN,N_PAR-skip);	
	/*INT_PREC=1.0e-5;
	save_res(((std::string)argv[1])+"/result.txt",min);	
	printf("***************************\n");
	printf("*** Second: eps=%.1e  ***\n",INT_PREC);
	printf("***************************\n");
	//ROOT::Minuit2::MnMigrad migrad2(theFCN, min.UserParameters() ,1);
	
	for(int i=0;i<20;i++){
		min=migrad(50,1);
		std::cout<<"Parameters "<<min.UserState()<<std::endl;
		std::cout<<std::scientific<<"fcn= "<<min.UserState().Fval()<<", edm= "<<min.UserState().Edm()<<"  "<<min.IsValid()<<std::endl;
		if(min.IsValid()){
			printf("FINE: %.3e/%.3e = %.3e\n", min.UserState().Edm(),  (min.UserState().Fval()),min.UserState().Edm()/ (min.UserState().Fval()));
			break;
		}
	}
	*/
	std::cout<<"Parameters "<<min.UserState()<<std::endl;
	std::cout<<"min= "<<min<<std::endl;
	std::fstream file;
	std::chrono::duration<double> time=walltime.now()-start;
	
	std::cout<<time.count()<<" seconds"<<std::endl;
	save_res(((std::string)argv[1])+"/result.txt",&min,&theFCN,N_PAR-skip);	
	
	
	
	return 0;
	
}
