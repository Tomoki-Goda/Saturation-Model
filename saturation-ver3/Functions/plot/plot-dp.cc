#include<iostream>
#include<fstream>
#include<cuba.h>
#include<cmath>
#include <vector>
#include <string>
#include"control.h"
#include"control-default.h"
#include"constants.h"
#include"Parameters.hh"

#include<chrono>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_sf.h>
#include<pthread.h>
#include<cuba.h>

extern double  INT_PREC;
#include"clenshaw.hh"
#include"gauss.hh"

#include"kt-formula.hh"
#include"options.h"


int N_APPROX=N_CHEB_R;

double INT_PREC=1.0e-4;

int main(int argc , char** argv){
	options opt=read_options(argc, argv);
	std::vector<double> param(10,0);

	char infilenames[500];
	char filenames[500];
	if(opt.path==""){
		opt.path=getenv("DIR");
	}
	if(opt.input_file_name!="result.txt"){	
		sprintf(infilenames,"%s",opt.input_file_name.c_str());
	}else{
		sprintf(infilenames,"%s/%s",opt.path.c_str(),opt.input_file_name.c_str());
	}
	
	printf("%s\n",infilenames);
	FILE* infile=fopen(infilenames,"r");
	read_parameters(infile,param);
	fclose(infile);
	double sigpar[10]={0},sudpar[10]={0};
	parameter(param,sigpar,sudpar);
	const double kt2max=1.0e+6;
////////////////////////////////////////////////////////////////////////////	
	#if GLUON_APPROX==1
		#if SIGMA_APPROX==-1||SIGMA_APPROX==0//if sigma is not in 1d grid
			SIGMA sigma;
			sigma.init(sigpar);
			INTEG dsigma(sigma);
			dsigma.init(sigpar,'l');
		#endif
		
		GLUON dipole_gluon(dsigma);
		dipole_gluon.init(sigpar);
	#else //only GBW K
		GLUON dipole_gluon;
		dipole_gluon.init(sigpar);
	#endif//GLUON_APPROX==1	
		
		
		
	 	printf("Initialized\n");
	 	//printf("Export\n");	
        if(opt.output_file_name!="out.txt"){
            sprintf(filenames,"%s",opt.output_file_name.c_str());
        }else{
		    sprintf(filenames,"%s/dp.txt",opt.path.c_str());
		}
		printf("%s\n",infilenames);
  		//FILE *outfile=fopen(filenames,"w");
//#if SUDAKOV>=1
		FILE *outfile=fopen(filenames,"w");
		//fclose(outfile);
		double val=0,x=0;
		double k2,mu2;
		const int x_len=100,y_len=100,y2_len=20;
		double arr[y2_len*y_len]={0};
		for(int k=0; k<x_len;++k){
			//x=X_MIN*pow(X_MAX/X_MIN,((double)k)/(x_len-1));
			x=1.0e-8*pow(1.0e+7,((double)k)/(x_len-1));
#if GLUON_APPROX==1
			dipole_gluon.set_x(x);
#endif
			printf("x= %.2e\n",x);
			
#pragma omp parallel private(k2,mu2)
{
			//double k2,mu2=0;
#pragma omp for schedule(dynamic)
			for(int j=0;j<y_len;++j){
				//k2=KT2_MIN*pow(kt2max/KT2_MIN,((double)j)/(y_len-1));
				k2=1.0e-2*pow(1.0e+3/1.0e-2,((double)j)/(y_len-1));
				printf("kt2= %.2e\n",k2);
#if SUDAKOV>=1
				for(int i=0;i<y2_len;++i){
					mu2=1.0e-1*pow(1.0e+5/1.0e-1,((double)i)/(y2_len-1));
					arr[j*y2_len+i]=dipole_gluon(x,k2,mu2);
				}
#else
				mu2=0;
				arr[j]=dipole_gluon(x,k2,mu2);
#endif
				printf("\033[1A\033[2K\r");

			}
}
			for(int j=0;j<y_len;++j){
				k2=1.0e-2*pow(1.0e+3/1.0e-2,((double)j)/(y_len-1));
/*#if SUDAKOV>=1
				for(int i=0;i<70;++i){
					mu2=1.0e-1*pow(1.0e+5/1.0e-1,((double)i)/69);
					val=arr[j*70+i];
					fprintf(outfile ,"%.10e\t%.10e\t%.10e\t%.10e\n",log(x),log(k2),log(mu2), val );
					//fprintf(outfile ,"%.10e\t%.10e\t%.10e\n",log(x),log(k2), val );
				}
#else*/	
				val=arr[j];
				fprintf(outfile ,"%.10e\t%.10e\t%.10e\n",log(x),log(k2),val/sigpar[0] );
//#endif
				

			}
			fprintf(outfile ,"\n");
			printf("\033[1A\033[2K\r");
			//gluon.set_max(kt2max,mu2);	
			
			//outfile=fopen(filenames,"a");
			//gluon.export_grid(outfile);
			//fclose(outfile);
		}
		fclose(outfile);
		printf("Gluon end\n");
//#else
//		gluon.set_max(kt2max);
//  		FILE *outfile=fopen(filenames,"w");
//		gluon.export_grid(outfile);
//		printf("Gluon end\n");
//	 	fclose(outfile);
//#endif		

	return 0;
}
