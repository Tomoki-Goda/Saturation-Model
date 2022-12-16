#include"Minuit2/FunctionMinimum.h"
#include"Minuit2/MnUserParameterState.h"
#include"Minuit2/MnPrint.h"
#include"Minuit2/MnMigrad.h"
#include"Minuit2/MnSimplex.h"
#include"Minuit2/MnHesse.h"
#include"Minuit2/FCNBase.h"
#include"Minuit2/MnMachinePrecision.h"


//#include<pthread.h>
#include"/home/tomoki/Numerics/clenshaw-curtis-gauss-legendre.hh"
#include<ctime>
//#include"./kt-formula.hh"
#ifndef R_FORMULA
	#define RFORMULA 0
#endif
//using namespace ROOT::Minuit2;

//namespace ROOT {
//namespace Minuit2 {
extern double F2_kt(double,double,double,double*);
extern double F2_r(double,double,double,double*);
class KtFCN : public ROOT::Minuit2::FCNBase {
	
	
	private:	
		//pthread_t threads[NF-2];
		//F2 F2_L(),F2_C(),F2_B();
		
		
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
		//virtual double operator()(const std::vector<double>&) const;
		//pthread_create(&(threads[k+j]), NULL,cckt,(void*)(&arg[j]));
		//double **ARGS;
		/*double F2(void* arg)const {
			double **args=(double**)arg;
			double x,Q2,mf2,*param;
			x=args[0][0];
			Q2=args[0][1];
			mf2=args[0][2];
			param=args[1];
			double res =F2_kt(x,Q2,mf2,param);
			*(args[3])=res;
			return(res);
		}*/
		
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
			static int licznik;
			licznik++;
			clock_t time;
			double x,Q2;
			double val;
			double chisq=0;
			int len=sizeof(par)/sizeof(par[0]);
			//double param[len];
			for(int i =0;i<len;i++){
				printf("%.5e\t",par[i]);
			}printf(" %d\n",licznik);
			//getchar();
			
			double sigpar[10],sudpar[10];
			parameter(par,sigpar, sudpar);
			
			//F2 F2L(MASS_L2, sigpar), F2C(MASS_C2, sigpar), F2B(MASS_B2, sigpar);
			//double *argsL[3],*argsC[3],*argsB[3];
			//double ext[3];
			//argsL[0]=ext;
			//argsC[0]=ext;
			//argsB[0]=ext;
			//double resL,resC,resB;
			//argsL[3]=&resL;
			//argsC[3]=&resC;
			//argsB[3]=&resB;
			time=clock();
			for(int i=0;i<MAX_N;i++){
				val=0;
				//ext[0]=X_DATA[i];
				//ext[1]=Q2_DATA[i];
				x=X_DATA[i];
				Q2=Q2_DATA[i];

				
#if R_FORMULA==1
				val+=(2.0/3.0)*F2_r( x,Q2,0,sigpar);
				//val+=F2_r( x,Q2,1,sigpar);
				val+=(4.0/9.0)*F2_r( x,Q2,2,sigpar);
				val+=(1.0/9.0)*F2_r( x,Q2,3,sigpar);
				//F2L((void*)argsL);
				//F2C((void*)argsC);
				//F2B((void*)argsB);
				//val=resL+resC+resB;
#else
				val+=(2.0/3.0)* F2_kt( x,Q2,MASS_L2,sigpar);
				val+=(4.0/9.0)* F2_kt( x,Q2,MASS_C2,sigpar);
				val+=(1.0/9.0)* F2_kt( x,Q2,MASS_B2,sigpar);		
#endif
				
				
				//printf("%.5e\t %.5e %d/%d in %f sec\n",val,(val-CS_DATA[i])/ERR_DATA[i],i,MAX_N,-((double)time)/CLOCKS_PER_SEC);
				
				chisq+=pow((val-CS_DATA[i])/ERR_DATA[i],2);
			}
			time-=clock();
			printf("CHISQ = %.5e (%.2f)  in %.1e sec\n",chisq, chisq/(MAX_N-len) , -((double)time)/CLOCKS_PER_SEC);			
			return chisq;
		}
		
		//void F2(void* par){
		//	double**param=(double**)par;
		//	
		 //return()
		//}

};
//}
//}
