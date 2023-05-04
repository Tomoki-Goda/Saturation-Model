#include"main.hh"

double double_round(double val,int i){
	int l=lrint(floor(log10(val)));
	l-=i-1;
	return(round( val*pow(10,-l) )*pow(10,l));
}
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
		 file<<min->UserState().Name(i)<<"\t"<<min->UserState().Value(i)<<"\t"<<min->UserState().Error(i)<<std::endl;
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


//////////////////////////////////////////////////////////
//
// MAIN 
//
////////////////////////////////////////////////////////////
int main(int argc, char** argv){
	
#if SCATTER==1
	std::fstream scatterfile; 
	scatterfile.open("home/tomoki/Saturation-Model/saturation-ver3/scatter.txt",std::fstream::out);
	scatterfile.close();
#endif
	std::chrono::system_clock walltime;
	std::chrono::time_point start= walltime.now();
	std::cout<<std::scientific<<std::endl;
	
	printf("********************Program Started********************.\n");
	printf("MODEL     = %d,               N_CHEB_R     = %d\n",MODEL,N_CHEB_R);
	printf("R_FORMULA = %d, MU02         = %d\n",R_FORMULA, MU02);
	printf("MODX      = %d, R_CHANGE_VAR = %d,             \t GLUON_APPROX = %d\n",MODX,R_CHANGE_VAR,GLUON_APPROX);
	printf("ALPHA_RUN = %d, \n", ALPHA_RUN);
	printf("FREEZE_QS2= %d, ADD_END      = %d,             \t THRESHOLD    = %d\n",FREEZE_QS2,ADD_END,THRESHOLD);
	printf("IBP          = %d\n",IBP);
	printf("ADJOINT   = %d, WW           = %d\n",ADJOINT,WW);
	printf("Directory = %s ,\t R  = [%.1e, %.1e]\n",(char*)argv[1],R_MIN, R_MAX);
	printf("SIGMA_APPROX = %d   N_CHEB   = %d              \t CHEB_D       = %d \n",SIGMA_APPROX,N_CHEB,CHEB_D);
	printf("*******************************************************.\n");
	

	KtFCN theFCN("/home/tomoki/Saturation-Model/saturation-ver3/data/hera_tot.dat",argv[1]);
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
////////////////////////////////////////////
// Set up the parameters
////////////////////////////////////////////
	for(unsigned i=0;i<N_PAR;i++){
#if (ALPHA_RUN==0||MU02!=0)
		if(par_name[i]=="mu02"){
			skip++;
			//std::cout<<"Skip "<< par_name[i]<<" = "<< MU02<<std::endl;;
			continue;
		}
#endif

		std::cout<<par_name[i]<<" = "<<par_start[i]<<std::endl;
		//std::cout<<N_PAR<<"  "<<skip <<" : parameter position="<<i<<std::endl;
		upar.Add(par_name[i], par_start[i],par_error[i]);
		upar.SetLimits(par_name[i],par_min[i],par_max[i]);//use migrad.removeLimits(<name>);
	}
#else//USE_RESULT==1
		printf("USING PREVIOUS RESULT %d\n",USE_RESULT);
	char resfile[100];
	sprintf(resfile,"%s/result.txt",argv[1]);
	
	char name[20];
	double ival,ierr;
	FILE* resinputfile=fopen(resfile,"r");
	fscanf(resinputfile,"%s %le",name,&ival);
	for(unsigned i=0;i<N_PAR;i++){
		fscanf(resinputfile,"%s %le %le",name,&ival,&ierr);
		if(strcmp(name,"chisq")==0){
			skip=N_PAR-i;
			break;
		}
#if USE_RESULT>0
		ival=double_round(ival,USE_RESULT);
#endif
		printf("%s %le %le \n",name,ival,ierr*10);
		upar.Add(name, ival, ival*5*pow(10,-USE_RESULT));
	}printf("\n");
	fclose(resinputfile);
#endif//USE_RESULT

//////////////////////////////////////////////////////////
// Start fitting
//////////////////////////////////////////////////////////
	ROOT::Minuit2::MnMachinePrecision prec;
	
	INT_PREC=5.0e-4 ;
	N_APPROX=N_CHEB_R/2;
	//prec.SetPrecision(INT_PREC);
	int flag=0;
	double goal=1;
	ROOT::Minuit2::MnSimplex simplex1(theFCN,upar,0);
	std::cout<<"TEST RUN 5, eps = "<<INT_PREC<<" N_APROX ="<<N_APPROX<<std::endl;	
	ROOT::Minuit2::FunctionMinimum min=simplex1(5,1);//Just initialization /check.
	std::cout<<"Parameters "<<min.UserState()<<std::endl;
	prec.SetPrecision(1.0e-6);
	
	INT_PREC=1.0e-3;
	N_APPROX=N_CHEB_R/3;
	//min=simplex1(25,1);
	for(int i=0;i<2;++i){
		printf("*****************************\n");
		printf("*** Simplex: eps=%.1e  N_APPROX=%d***\n",(double)INT_PREC,N_APPROX);
		printf("*****************************\n");
		//min=simplex1(100,pow(10,1-i));
//		if(i==1){		
			for(int k=0;k<N_PAR-skip;++k){
				//min=simplex1(50,0);
				simplex1.Fix(k);
//#if MODEL==1
//				for(int j=0;j<N_PAR-skip;++j){
//					if(i==j){
//						continue;
//					}
//					simplex1.Fix(j);
//					min=simplex1(100,pow(10,1-i));
//					std::cout<<"Parameters "<<min.UserState()<<std::endl;
//					simplex1.Release(j);
//				}
//#endif
				min=simplex1(50,50/(k+1));
				std::cout<<"Parameters "<<min.UserState()<<std::endl;
				simplex1.Release(k);
			}
//		}
		
		min=simplex1(50,25);
		INT_PREC/=2;
		N_APPROX=(int)(N_APPROX*2);
		save_res(((std::string)argv[1])+"/result.txt",&min,&theFCN,N_PAR-skip);
	}

/////////////////////////////////////////////////////////////////
// Try MiGrad
//////////////////////////////////////////////////////////////////	
	std::vector<int> fix;
	std::vector<int> free;
	fix.reserve(10);
	free.reserve(10);
	
	INT_PREC=2.5e-4;
	N_APPROX=(2*N_CHEB_R)/3;
	//prec.SetPrecision(INT_PREC);
	printf("***************************\n");
	printf("*** First: eps=%.1e  N_APROX=%d***\n",(double)INT_PREC,N_APPROX);
	printf("***************************\n");
{	
	ROOT::Minuit2::MnMigrad migrad(theFCN,min.UserParameters(),0);
	for(int i=0;i<(N_PAR-skip);i++ ){
		migrad.RemoveLimits(i);
	}
	min=migrad(20,25);
}
	
	
	std::cout<<"Parameters "<<min.UserState()<<std::endl;
	goal=25;
	for(int k=0;k<5;k++){
		ROOT::Minuit2::MnMigrad migrad(theFCN,min.UserParameters(),0);
		//////// try fixing////////////
		fix.clear();
		free.clear();
		for(int ii=0;ii<N_PAR-skip;++ii){
			free.push_back(ii);
		}
		for(int j=0;j<N_PAR-skip-1;++j){
			std::vector<double> cor=min.UserState().GlobalCC().GlobalCC();
			for(int ii=0;ii<cor.size();++ii){
				printf("cor %d = %f\n",ii,cor[ii]);
			}
			int fixone=std::distance(cor.begin(),std::max_element(cor.begin(),cor.end()));
			
			printf("Fix %d (st/nd/th) pos (%d) \n",free[fixone],fixone);
			
			migrad.Fix(free[fixone]);
			
			fix.push_back(free[fixone]);
			free.erase(free.begin()+fixone);
			
			min=migrad(20,goal);
			
			std::cout<<"Parameters "<<min.UserState()<<std::endl;
			printf("Cov= %d\n",min.UserState().CovarianceStatus() );
			printf("Valid: %d \tCovariance: %d\n",min.IsValid(),min.HasCovariance());
			if(min.IsValid()){
				printf("Fixing successful\n");
				break;
			}else{
				printf("Fixing insufficient\n");
			}
		}
		
		for(int ii=0;ii<fix.size();++ii){
			migrad.Release(fix[ii]);
		}
		
//////////////////////////////////////////////////////////////////////	
//		for(int i=0;i<N_PAR-skip;++i){
			min=migrad(50,goal);
			std::cout<<"Parameters "<<min.UserState()<<std::endl;
//			if(min.IsValid()){
//				break;
//			}
//		}
		
		
		printf("Cov= %d\n",min.UserState().CovarianceStatus() );
		printf("Valid: %d \tCovariance: %d\n",min.IsValid(),min.HasCovariance());
		
		save_res(((std::string)argv[1])+"/result.txt",&min,&theFCN,N_PAR-skip);	
		if(min.IsValid()){
			printf(" %.3e/%.3e = %.3e\n", min.UserState().Edm(),  (min.UserState().Fval()),min.UserState().Edm()/ (min.UserState().Fval()));
			break;
		}
		INT_PREC/=1.4142;//sqrt(2)
		N_APPROX=((int)(1.4142*N_APPROX));

		std::cout<<"Parameters "<<min.UserState()<<std::endl;
		
		
		//////////Try with simplex///////////////
		//ROOT::Minuit2::MnSimplex simplex(theFCN,min.UserParameters(),0);
		//min=simplex(25,goal);
		//std::cout<<"Parameters "<<min.UserState()<<std::endl;
	}
	
	INT_PREC=1.0e-4;
	N_APPROX=N_CHEB_R;

	printf("***************************\n");
	printf("*** Second: eps=%.1e  N_APPROX=%d***\n",(double)INT_PREC,N_APPROX);
	printf("***************************\n");
	goal=10;
	for(int k=0;k<5;k++){
		ROOT::Minuit2::MnMigrad migrad(theFCN,min.UserParameters(),1);
		//////// try fixing////////////
		fix.clear();
		free.clear();
		for(int ii=0;ii<N_PAR-skip;++ii){
			free.push_back(ii);
		}
		for(int j=0;j<N_PAR-skip-1;++j){
			std::vector<double> cor=min.UserState().GlobalCC().GlobalCC();
			for(int ii=0;ii<cor.size();++ii){
				printf("cor %d = %f\n",ii,cor[ii]);
			}
			int fixone=std::distance(cor.begin(),std::max_element(cor.begin(),cor.end()));
			
			printf("Fix %d (st/nd/th) pos (%d) \n",free[fixone],fixone);
			
			migrad.Fix(free[fixone]);
			
			fix.push_back(free[fixone]);
			free.erase(free.begin()+fixone);
			
			min=migrad(20,goal);
			
			std::cout<<"Parameters "<<min.UserState()<<std::endl;
			printf("Cov= %d\n",min.UserState().CovarianceStatus() );
			printf("Valid: %d \tCovariance: %d\n",min.IsValid(),min.HasCovariance());
			if(min.IsValid()){
				printf("Fixing successful\n");
				break;
			}else{
				printf("Fixing insufficient\n");
			}
		}
		
		for(int ii=0;ii<fix.size();++ii){
			migrad.Release(fix[ii]);
		}
//		for(int i=0;i<N_PAR-skip;++i){
			min=migrad(50,goal);
			std::cout<<"Parameters "<<min.UserState()<<std::endl;
//		}

		if(min.IsValid()&&(min.UserState().CovarianceStatus()==3 )){
			printf(" %.3e/%.3e = %.3e\n", min.UserState().Edm(),  (min.UserState().Fval()),min.UserState().Edm()/ (min.UserState().Fval()));
			break;
		}

		INT_PREC/=1.41421356;//sqrt(2)
		N_APPROX=((int)(1.41421356*N_APPROX));
				
		std::cout<<"Parameters "<<min.UserState()<<std::endl;
		//
		//ROOT::Minuit2::MnSimplex simplex(theFCN,min.UserParameters(),0);
		//min=simplex(30,goal);	
	}

	
	

	std::cout<<"Parameters "<<min.UserState()<<std::endl;
	std::cout<<"min= "<<min<<std::endl;
	std::chrono::duration<double> time=walltime.now()-start;
	std::cout<<time.count()<<" seconds"<<std::endl;
	return 0;
	
}

	
	
	
	
