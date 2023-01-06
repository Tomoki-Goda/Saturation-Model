#include"Minuit2/FunctionMinimum.h"
#include"Minuit2/MnUserParameterState.h"
#include"Minuit2/MnPrint.h"
#include"Minuit2/MnMigrad.h"
#include"Minuit2/MnSimplex.h"
#include"Minuit2/MnHesse.h"
#include"Minuit2/MnEigen.h"
#include"Minuit2/FCNBase.h"
#include"Minuit2/MnMachinePrecision.h"

//#include"/home/tomoki/Numerics/clenshaw-curtis-gauss-legendre.hh"
#include"./clenshaw-curtis.hh"
#include<ctime>
#include<chrono>

#ifndef R_FORMULA
	#define RFORMULA 0
#endif

extern double F2_kt(const  PREC,const  PREC,const PREC, const PREC(&)[]);
//extern double F2_r(const double,double,double,double*);
class KtFCN : public ROOT::Minuit2::FCNBase {
	
	
	private:	
		double *__restrict Q2_DATA=(double*)malloc(0);
		double *__restrict X_DATA=(double*)malloc(0);
		double *__restrict CS_DATA=(double*)malloc(0);
		double *__restrict ERR_DATA=(double*)malloc(0);
		
		void data_alloc(int n){
			MAX_N=n;
			Q2_DATA=(double*  )realloc(Q2_DATA,MAX_N*sizeof(double));
			X_DATA=(double* )realloc(X_DATA,MAX_N*sizeof(double));
			CS_DATA=(double*)realloc(CS_DATA,MAX_N*sizeof(double));
			ERR_DATA=(double*)realloc(ERR_DATA,MAX_N*sizeof(double));
		}
		double Up() const {return 1;}
		
	public:
		unsigned MAX_N=0;
		explicit KtFCN(std::string data_file ){
		//KtFCN(char* data_file){
			data_alloc(500);
			const double alpha =1.0/137 ;//fine structure const 1/137;
			const double xmp0 = 0.93827;//proton mass in GeV
			const double units =1.0/389.40; //2.56819e-3; //micro-barn to GeV^-2
			unsigned i=0;
			unsigned j=0;
			double fac;
			//std::fstream file;
			//file.open(data_file, std::fstream::in);
			FILE* file=fopen(data_file.c_str(),"r");
			double wdata;
			//while((!file.eof())&&(j<597)){
			if(file==NULL){
				printf("file not found.\n");
			}
			while((!feof(file))&&(j<597)){
				fscanf(file,"%lf %lf %lf %lf %lf", (Q2_DATA+i),(X_DATA+i),&wdata,(CS_DATA+i),(ERR_DATA+i)); 
				//printf("%le %le %le %le %le\n", *(Q2_DATA+i),*(X_DATA+i),wdata,*(CS_DATA+i),*(ERR_DATA+i)); 
				
			
				if((X_DATA[i]<=X_MAX)&&( Q2_DATA[i]<=Q2_MAX)){
					//printf("%le %le %le %le %le\n", *(Q2_DATA+i),*(X_DATA+i),wdata,*(CS_DATA+i),*(ERR_DATA+i)); 
					/////formula in I. Abt et al 2017////////
					fac = pow(Q2_DATA[i],2)*(1-X_DATA[i])/ (4*pow(PI,2)*alpha*(Q2_DATA[i]+pow(2*X_DATA[i]*xmp0,2)));
					fac*=units;	
					CS_DATA[i] = fac*CS_DATA[i];
					ERR_DATA[i] = fac*ERR_DATA[i];
					//fprintf(stdout, "%d: %lE %lE %lE %lE %lE\n",i+1,*(X_DATA+i), *(Y_DATA+i), *(Q2_DATA+i), *(CS_DATA+i), *(ERR_DATA+i));
					//printf("%.3e %.3e %.3e %.3e %.3e\n", Q2_DATA[i],X_DATA[i],wdata,CS_DATA[i],ERR_DATA[i]); 
					++i;
				}
				if(i==MAX_N){
					data_alloc(MAX_N*2);
				}	
				
				++j;
			}
			MAX_N=i;	
			//file.close();	
			fclose(file);	
			printf("HERA data loaded\n");
		
		}
		~KtFCN(){
			free(Q2_DATA);
			free(X_DATA);
			free(CS_DATA);
			free(ERR_DATA);
		}
		
		double operator()(const std::vector<double>& par)const{
			std::chrono::system_clock walltime;
			std::chrono::time_point start= walltime.now();

			clock_t time=clock();
			static int licznik;
			++licznik;
			double x,Q2;
			double val;
			double chisq=0;
			const int len=sizeof(par)/sizeof(par[0]);
			//double param[len];
#if PRINT_PROGRESS!=0
			if((licznik/PRINT_PROGRESS)*PRINT_PROGRESS==licznik){			
				for(int i =0;i<len;++i){
					printf("%.5e\t",par[i]);
				}printf(" %d\n",licznik);
			//getchar();
			}
#endif			
			PREC sigpar[10],sudpar[10];
			parameter(par,sigpar, sudpar);
			for(int i=0;i<MAX_N;++i){
				val=0;
				//x=X_DATA[i];
				//Q2=Q2_DATA[i];
				val=F2_kt(X_DATA[i],Q2_DATA[i],0,sigpar);//summation over flavour is done at the level of integrand.
				
				chisq+=pow((val-CS_DATA[i])/ERR_DATA[i],2);
			}

			std::chrono::duration<double> interval=walltime.now()-start;
			time-=clock();
#if PRINT_PROGRESS!=0
			if((licznik/PRINT_PROGRESS)*PRINT_PROGRESS==licznik){
				printf("CHISQ = %.5e (%.2f) \t",chisq, chisq/(MAX_N-len) );//, -((double)time)/CLOCKS_PER_SEC);			
				std::cout<<interval.count()<<" seconds, ("<< -((double)time)/CLOCKS_PER_SEC<<" CPU seconds)"<<std::endl;
			}
#endif	
			if(std::isnan(chisq)+std::isinf(chisq)){
				printf("theFCN:: %.3e encountered\n", chisq);
				return(0);
			}
			return ((double)chisq);
		}
		

};

ROOT::Minuit2::FunctionMinimum migrad(const ROOT::Minuit2::FCNBase& theFCN, ROOT::Minuit2::MnUserParameterState& stat, const int  str){
		ROOT::Minuit2::MnStrategy strat(str);
		ROOT::Minuit2::MnMigrad migrad2(theFCN,stat,strat);	
		ROOT::Minuit2::FunctionMinimum min=migrad2(50,10);
		ROOT::Minuit2::MnEigen eigen;

		std::cout <<"Eigen values: " ;
		for (double val: eigen(min.UserCovariance())){
   			std::cout << val << ' ';
   		}
   		std::cout<<"\n"<<std::endl;
   		std::cout<<"Parameters "<<min.UserState()<<"\n"<<std::endl; 
   		stat=min.UserState();
		return min;		
}
	
	

