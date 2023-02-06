#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<vector>
#include"./control.h"
#include"./control-default.h"
#include"./constants.h"
#include"./Parameters.hh"
#include"cfortran.h"

#ifndef USE_RESULT
	#define USE_RESULT 0
#endif


PREC INT_PREC=DGAUSS_PREC;
int N_APPROX=N_CHEB_R;
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
	file<<"Cov\t"<<min->UserState().CovarianceStatus()<<std::endl;
	file<<"eps\t"<<INT_PREC<<std::endl;
	
	file.close();
	return 0;
}

int main(int argc, char** argv){
	
#if SCATTER==1
	std::fstream scatterfile; 
	scatterfile.open("home/tomoki/Saturation-Model/saturation-ver3/scatter.txt",std::fstream::out);
	scatterfile.close();
#endif
	std::chrono::system_clock walltime;
	std::chrono::time_point start= walltime.now();
	std::cout<<std::scientific<<std::endl;
	
	printf("Program Started.\n");
	printf("Directory=%s\n",(char*)argv[1]);
	

	KtFCN theFCN("/home/tomoki/Saturation-Model/saturation-ver3/data/hera_tot.dat");
	ROOT::Minuit2::MnUserParameters upar;

#if MU202!=0
#if (MODEL==3 && INDEPENDENT_RMAX==0)
	printf(" mu202 is not free. Cannot be controlled\n");
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
#if USE_RESULT==0
//	for(unsigned i=0;(i-skip)<N_PAR;i++){
	for(unsigned i=0;i<N_PAR;i++){
#if (ALPHA_RUN==0||MU02!=0)
		if(par_name[i]=="mu102"){
			skip++;
			//std::cout<<"Skip "<< par_name[i]<<" = "<< MU02<<std::endl;;
			continue;
		}
#endif
		std::cout<<par_name[i]<<" = "<<par_start[i]<<std::endl;
		std::cout<<N_PAR<<"  "<<skip <<" : parameter position="<<i<<std::endl;
		upar.Add(par_name[i], par_start[i],par_error[i]);
		upar.SetLimits(par_name[i],par_min[i],par_max[i]);//use migrad.removeLimits(<name>);
	}
#else//USE_RESULT==1
	char resfile[100];
	sprintf(resfile,"%s/result.txt",argv[1]);
	
	char name[20];
	double ival,ierr;
	FILE* resinputfile=fopen(resfile,"r");
	fscanf(resinputfile,"%s %le",name,&ival);
	for(unsigned i=0;i<N_PAR;i++){
		fscanf(resinputfile,"%s %le %le",name,&ival,&ierr);
		printf("%s %le %le \n",name,ival,ierr);
		if(strcmp(name,"chisq")==0){
			skip=N_PAR-i;
			break;
		}
		upar.Add(name, ival,ierr);
	}printf("\n");
	fclose(resinputfile);
#endif//USE_RESULT
	ROOT::Minuit2::MnMachinePrecision prec;
	//prec.SetPrecision(1.0e-8);
	INT_PREC=5.0e-4;
	N_APPROX=N_CHEB_R/2;
	//prec.SetPrecision(INT_PREC);
	int flag=0;
	double goal=1;
	ROOT::Minuit2::MnSimplex simplex1(theFCN,upar,0);
	std::cout<<"TEST RUN 15, eps = "<<INT_PREC<<" N_APROX ="<<N_APPROX<<std::endl;	
	ROOT::Minuit2::FunctionMinimum min=simplex1(15,1);//Just initialization /check.
	ROOT::Minuit2::FunctionMinimum min_prev=min;
	//ROOT::Minuit2::MnEigen eigen;
	//min_prev=min;
	//std::cout<<"Parameters "<<min_prev.UserState()<<std::endl;
	std::cout<<"Parameters "<<min.UserState()<<std::endl;
	INT_PREC=1.0e-3;
	N_APPROX=N_CHEB_R/4;
	//ROOT::Minuit2::MnSimplex simplex2(theFCN,upar,0);
	
	for(int i=0;i<2;++i){
		prec.SetPrecision(INT_PREC*2);
		printf("*****************************\n");
		printf("*** Simplex: eps=%.1e  N_APPROX=%d***\n",(double)INT_PREC,N_APPROX);
		printf("*****************************\n");
		min=simplex1(100,1);
		INT_PREC/=2;
		//N_APPROX*=2;
		save_res(((std::string)argv[1])+"/result.txt",&min,&theFCN,N_PAR-skip);
	}
	
	
//	std::cout<<upar<<std::endl;
//	printf("N_PAR-skip=%d-%d=%d\n",N_PAR,skip,N_PAR-skip );
	INT_PREC=5.0e-4;
	N_APPROX=N_CHEB_R/2;
	prec.SetPrecision(INT_PREC*2);
	printf("***************************\n");
	printf("*** First: eps=%.1e  N_APROX=%d***\n",(double)INT_PREC,N_APPROX);
	//printf("N_PAR-skip=%d\n",N_PAR-skip );
	//std::cout<<upar<<std::endl;
	printf("***************************\n");
	ROOT::Minuit2::MnHesse hesse;
	ROOT::Minuit2::MnUserParameterState stat=hesse(theFCN,min.UserParameters());
//	ROOT::Minuit2::MnUserParameterState stat=hesse(theFCN,upar);
	for(int i=0;i<(N_PAR-skip);i++ ){
		stat.RemoveLimits(i);
	}
	
	std::cout<<"Hesse "<<stat<<std::endl; 
	ROOT::Minuit2::MnStrategy strat0(0);  	
	ROOT::Minuit2::MnMigrad migrad1(theFCN,stat,strat0);
	
//	ROOT::Minuit2::FunctionMinimum	min=migrad1(10,goal);
   	goal=10;
	for(int i=0;i<5;i++){
		min=migrad1(10*(i+1),goal);
				
		printf("EDM/FVal %.3e/%.3e = %.3e\n",min.UserState().Edm(),min.UserState().Fval(),((min.UserState().Edm())/(min.UserState().Fval())) );
		printf("Cov= %d\n",min.UserState().CovarianceStatus() );
		printf("Valid: %d \tCovariance: %d\n",min.IsValid(),min.HasCovariance());
		
		save_res(((std::string)argv[1])+"/result.txt",&min,&theFCN,N_PAR-skip);	
		if(min.IsValid()&&(min.UserState().CovarianceStatus()==3 )){
			printf(" %.3e/%.3e = %.3e\n", min.UserState().Edm(),  (min.UserState().Fval()),min.UserState().Edm()/ (min.UserState().Fval()));
			break;
		}
		
	}
	
	//save_res(((std::string)argv[1])+"/result.txt",&min,&theFCN,N_PAR-skip);	

	INT_PREC=1.0e-4;
	N_APPROX=N_CHEB_R;
	prec.SetPrecision(4*INT_PREC);
	printf("***************************\n");
	printf("*** Second: eps=%.1e  N_APPROX=%d***\n",(double)INT_PREC,N_APPROX);
	printf("***************************\n");
	stat=hesse(theFCN,min.UserParameters());
	std::cout<<"Hesse "<<stat<<std::endl;   	
	ROOT::Minuit2::MnStrategy strat1(1);  	
	ROOT::Minuit2::MnMigrad migrad2(theFCN,stat,strat1);
	
   	goal=5;
	for(int i=0;i<5;i++){
		min=migrad2(10*(i+1),goal);
				
		printf("EDM/FVal %.3e/%.3e = %.3e\n",min.UserState().Edm(),min.UserState().Fval(),((min.UserState().Edm())/(min.UserState().Fval())) );
		printf("Cov= %d\n",min.UserState().CovarianceStatus() );
		printf("Valid: %d \tCovariance: %d\n",min.IsValid(),min.HasCovariance());
		
		save_res(((std::string)argv[1])+"/result.txt",&min,&theFCN,N_PAR-skip);
		if(min.IsValid()&&(min.UserState().CovarianceStatus()==3 )){
			printf(" %.3e/%.3e = %.3e\n", min.UserState().Edm(),  (min.UserState().Fval()),min.UserState().Edm()/ (min.UserState().Fval()));
			break;
		}
		
	}
	//save_res(((std::string)argv[1])+"/result.txt",&min,&theFCN,N_PAR-skip);	
	
	std::cout<<"Parameters "<<min.UserState()<<std::endl;
	std::cout<<"min= "<<min<<std::endl;
	//std::fstream file;
	std::chrono::duration<PREC> time=walltime.now()-start;
	
	std::cout<<time.count()<<" seconds"<<std::endl;
	//save_res(((std::string)argv[1])+"/result.txt",&min,&theFCN,N_PAR-skip);	
	return 0;
	
}

	
	
	
	
