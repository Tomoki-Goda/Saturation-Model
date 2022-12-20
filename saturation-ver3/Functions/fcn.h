#include"Minuit2/FunctionMinimum.h"
#include"Minuit2/MnUserParameterState.h"
#include"Minuit2/MnPrint.h"
#include"Minuit2/MnMigrad.h"
#include"Minuit2/MnSimplex.h"
#include"Minuit2/MnHesse.h"
#include"Minuit2/FCNBase.h"
#include"Minuit2/MnMachinePrecision.h"

#include"/home/tomoki/Numerics/clenshaw-curtis-gauss-legendre.hh"
#include<ctime>
#include<chrono>

#ifndef R_FORMULA
	#define RFORMULA 0
#endif

extern double F2_kt(const double,const double,const double, const double*);
//extern double F2_r(const double,double,double,double*);
class KtFCN : public ROOT::Minuit2::FCNBase {
	
	
	private:	
		double *Q2_DATA=(double*)malloc(0);
		double *X_DATA=(double*)malloc(0);
		double *CS_DATA=(double*)malloc(0);
		double *ERR_DATA=(double*)malloc(0);
		
		void data_alloc(int n){
			MAX_N=n;
			Q2_DATA=(double*)realloc(Q2_DATA,MAX_N*sizeof(double));
			X_DATA=(double*)realloc(X_DATA,MAX_N*sizeof(double));
			CS_DATA=(double*)realloc(CS_DATA,MAX_N*sizeof(double));
			ERR_DATA=(double*)realloc(ERR_DATA,MAX_N*sizeof(double));
		}
		double Up() const {return 1;}
		
	public:
		unsigned MAX_N=600;
		KtFCN(std::string data_file ){
		//KtFCN(char* data_file){
			data_alloc(500);
			double alpha =1.0/137 ;//fine structure const 1/137;
			double xmp0 = 0.93827;//proton mass in GeV
			double units =1.0/389.40; //2.56819e-3; //micro-barn to GeV^-2
			unsigned i=0;
			unsigned j=0;
			double fac;
			//std::fstream file;
			//file.open(data_file, std::fstream::in);
			FILE* file=fopen(data_file.c_str(),"r");
			double wdata;
			//while((!file.eof())&&(j<597)){
			while((!feof(file))&&(j<597)){
				fscanf(file,"%lE %lE %lE %lE %lE", (Q2_DATA+i),(X_DATA+i),&wdata,(CS_DATA+i),(ERR_DATA+i)); 
				/////formula in I. Abt et al 2017////////
				fac = pow(Q2_DATA[i],2)*(1-X_DATA[i])/ (4*pow(PI,2)*alpha*(Q2_DATA[i]+pow(2*X_DATA[i]*xmp0,2)));
				fac*=units;	
				CS_DATA[i] = fac*CS_DATA[i];
				ERR_DATA[i] = fac*ERR_DATA[i];
			
				if((X_DATA[i]<=X_MAX)&&( Q2_DATA[i]<=Q2_MAX)){
					//fprintf(stdout, "%d: %lE %lE %lE %lE %lE\n",i+1,*(X_DATA+i), *(Y_DATA+i), *(Q2_DATA+i), *(CS_DATA+i), *(ERR_DATA+i));
					i++;
				}
				if(i==MAX_N){
					data_alloc(MAX_N*2);
				}	
				
				j++;
			}
			MAX_N=i;	
			//file.close();	
			fclose(file);	
		
		}
		~KtFCN(){
			free(Q2_DATA);
			free(X_DATA);
			free(CS_DATA);
			free(ERR_DATA);
		}
		
		double operator()(const std::vector<double> & par)const{
			std::chrono::system_clock walltime;
			std::chrono::time_point start= walltime.now();

			clock_t time=clock();
			static int licznik;
			licznik++;
			double x,Q2;
			double val;
			double chisq=0;
			int len=sizeof(par)/sizeof(par[0]);
			//double param[len];
#if PRINT_PROGRESS!=0
			if((licznik/PRINT_PROGRESS)*PRINT_PROGRESS==licznik){			
				for(int i =0;i<len;i++){
					printf("%.5e\t",par[i]);
				}printf(" %d\n",licznik);
			//getchar();
			}
#endif			
			double sigpar[10],sudpar[10];
			parameter(par,sigpar, sudpar);
			for(int i=0;i<MAX_N;i++){
				val=0;
				x=X_DATA[i];
				Q2=Q2_DATA[i];
				val=F2_kt(x,Q2,0,sigpar);//summation over flavour is done at the level of integrand.
				
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
			return chisq;
		}
		

};
