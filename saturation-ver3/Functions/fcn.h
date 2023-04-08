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
//#include"./clenshaw-curtis.hh"
//#include"./kt-formula.hh"

//extern double F2_kt(const  double,const  double,const double, const double(&)[]);
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
		FILE* file=NULL;
		std::string directory;
	public:
		unsigned MAX_N=0;
		int flag=0;
		explicit KtFCN(std::string data_file,std::string dir ){
			
				//std::string file_name=dir+"/log.txt";
				//file=fopen(file_name.c_str(),"a");
				directory=dir;
			

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
				exit(1);
			}
			while((!feof(file))&&(j<597)){
				fscanf(file,"%lf %lf %lf %lf %lf", (Q2_DATA+i),(X_DATA+i),&wdata,(CS_DATA+i),(ERR_DATA+i)); 
				if((X_DATA[i]<=X_DATA_MAX)&&( Q2_DATA[i]<=Q2_MAX)){
					fac = pow(Q2_DATA[i],2)*(1-X_DATA[i])/ (4*pow(PI,2)*alpha*(Q2_DATA[i]+pow(2*X_DATA[i]*xmp0,2)));
					fac*=units;	
					CS_DATA[i] = fac*CS_DATA[i];
					ERR_DATA[i] = fac*ERR_DATA[i];
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
			//fclose(file);
		}
		double operator()(const std::vector<double>& par)const{
			std::chrono::system_clock walltime;
			int cub=0;
			cubacores(&cub,&cub);
				
			std::chrono::time_point start= walltime.now();

			clock_t time=clock();
			static int licznik;
			++licznik;
			double x,Q2;
			
			double chisq=0;
			const int len=par.size();
			//printf("%d parameters\n",len);
			//double param[len];
#if PRINT_PROGRESS!=0
			if((licznik/PRINT_PROGRESS)*PRINT_PROGRESS==licznik){			
				for(int i =0;i<len;++i){
					printf("%.5e\t",par[i]);
				}printf(" %d\n",licznik);
			//getchar();
			}
#endif			
			double sigpar[10]={0},sudpar[10]={0};
			parameter(par,sigpar, sudpar);//Format
//////////////////////////////////////////////////////////////////////////////////////////
// Below is very messy
//////////////////////////////////////////////////////////////////////////////////////////


#if R_FORMULA==0///////////////////////////////////////////////////////////////////////
	#if GLUON_APPROX==1
		#if SIGMA_APPROX==-1||SIGMA_APPROX==0//if sigma is not in 1d grid
			SIGMA sigma;
			sigma.init(sigpar);
			DSIGMA dsigma(sigma);
			dsigma.init(sigpar,'l');
 		#else //-2 and 1  
 			SIGMA sigma;
			sigma.init(sigpar);
			DSIGMA dsigma(sigma);
			dsigma.init(N_APPROX+250,sigpar,'l');
		#endif
		
		GLUON dipole_gluon(dsigma);
		dipole_gluon.init(sigpar);
		Approx_aF<GLUON> gluon(dipole_gluon);
		gluon.init(N_APPROX+100,N_APPROX+100,sigpar);
		const double kt2max=9.0e+4;
		gluon.set_max(kt2max);
	#else               //only GBW K
		GLUON gluon;
		gluon.init(sigpar);
	#endif//GLUON_APPROX==1	
#endif//R_FORMULA
///////////////////////////////////////////////////////////////////////////////////////////		
			double arr[MAX_N];
			double *arr1,*arr2;
			if(flag==1){
				arr1=(double*)malloc(MAX_N*sizeof(double));
				arr2=(double*)malloc(MAX_N*sizeof(double));
			}			
#pragma omp parallel 
{ 

//////////////////////////////////////////////////////////////////////////////////////////
#if R_FORMULA==1	
	SIGMA sigma[3]={SIGMA() ,SIGMA() ,SIGMA() };
		sigma[0].init(sigpar);
		sigma[1].init(sigpar);
		sigma[2].init(sigpar);	

	#if SIGMA_APPROX==-2||SIGMA_APPROX==1 // AS GBW/BGK K   This require instance per thread
			DSIGMA dsigma[3]={DSIGMA(sigma[0]) ,DSIGMA(sigma[1]) ,DSIGMA(sigma[2]) };
			dsigma[0].init(N_APPROX+250,sigpar,'s');
			dsigma[1].init(N_APPROX+250,sigpar,'s');
			dsigma[2].init(N_APPROX+250,sigpar,'s');
			
			Integrand_r<DSIGMA> integrands[3]={
				Integrand_r<DSIGMA>(dsigma[0]) ,
				Integrand_r<DSIGMA>(dsigma[1]) ,
				Integrand_r<DSIGMA>(dsigma[2])
			};
			F2_kt<Integrand_r<DSIGMA>> F2(integrands);
	#elif SIGMA_APPROX==0||SIGMA_APPROX==-1
		SIGMA (&dsigma)[3]=sigma;
			Integrand_r<SIGMA> integrands[3]={
				Integrand_r<SIGMA>(dsigma[0]) ,
				Integrand_r<SIGMA>(dsigma[1]) ,
				Integrand_r<SIGMA>(dsigma[2])
			};
			F2_kt<Integrand_r<SIGMA>> F2(integrands);
	#endif///////////////////////////////////////////////////////////////////////////////////

#else//R_FORMULA==0
			Integrand_kt<Approx_aF<GLUON>> integrands[3]={
				Integrand_kt<Approx_aF<GLUON>>( gluon),
				Integrand_kt<Approx_aF<GLUON>>( gluon),
				Integrand_kt<Approx_aF<GLUON>>( gluon)
			};
			F2_kt<Integrand_kt<Approx_aF<GLUON>>> F2(integrands);
#endif//R_FORMULA

#pragma omp for schedule(dynamic)
			for(int i=0;i<MAX_N;++i){
				//F2_kt F2(sigpar);
				arr[i]=F2(X_DATA[i],Q2_DATA[i],0);//summation over flavour is done at the level of integrand.
				if(flag==1){
					arr1[i]=F2(X_DATA[i]*0.75,Q2_DATA[i],0);//don't forget to match fprintf below
					arr2[i]=F2(X_DATA[i]*1.25,Q2_DATA[i],0);
				}
				if(i>0){
					printf("\033[1A \033[2K");
				}
				printf("%d: val=%.2e data= %.2e chisq=%.2e x=%.2e Q2=%.2e\n",i,arr[i],CS_DATA[i],pow(fabs(arr[i]-CS_DATA[i])/ERR_DATA[i],2),X_DATA[i],Q2_DATA[i]);
			}
}
			chisq=0;
			for(int i=0;i<MAX_N;++i){
				chisq+=pow((arr[i]-CS_DATA[i])/ERR_DATA[i],2);
			}
			if(flag==1){	
				//FILE* file0=fopen((directory+"/aF.txt").c_str(),"w" );
				FILE* file1=fopen((directory+"/data.txt").c_str(),"w" );
				FILE* file2=fopen((directory+"/F2.txt").c_str(),"w" );
				//gluon.export_grid(file0);
				for(int i=0;i<MAX_N;++i){
					fprintf(file1, "%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n",X_DATA[i],Q2_DATA[i],CS_DATA[i],ERR_DATA[i],arr[i]);
					
					fprintf(file2, "%.5e\t%.5e\t%.5e\n",X_DATA[i]*0.75,Q2_DATA[i],arr1[i]);
					fprintf(file2, "%.5e\t%.5e\t%.5e\n",X_DATA[i],Q2_DATA[i],arr[i]);
					fprintf(file2, "%.5e\t%.5e\t%.5e\n",X_DATA[i]*1.25,Q2_DATA[i],arr2[i]);
				}
				//fclose(file0);
				fclose(file1);
				fclose(file2);
				free(arr1);
				free(arr2);
			}else{
				//if(file==NULL){
					FILE* file=fopen((directory+"/log.txt").c_str(),"a");
				//}
				for(int i =0;i<len;++i){
					fprintf(file,"%.5e\t",par[i]);
				}
				fprintf(file,"%.5e\n",chisq);
				fflush(file);
				fclose(file);
			}
			printf("\033[1A \033[2K");
			
			
			std::chrono::duration<double> interval=walltime.now()-start;
			time-=clock();
#if PRINT_PROGRESS!=0
			if((licznik/PRINT_PROGRESS)*PRINT_PROGRESS==licznik){
				printf("CHISQ = %.5e (%.3f) \t",chisq, chisq/(MAX_N-len) );//, -((double)time)/CLOCKS_PER_SEC);	
				printf("%d, %.2e\t",N_APPROX,INT_PREC );		
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


