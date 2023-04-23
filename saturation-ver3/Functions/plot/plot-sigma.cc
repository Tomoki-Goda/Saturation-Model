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

double INT_PREC=1.0e-5;


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
////////////////////////////////////////////////////////////////////////////	
			SIGMA sigma;
			sigma.init(sigpar);
        
		if(opt.output_file_name!="out.txt"){
            sprintf(filenames,"%s",opt.output_file_name.c_str());
        }else{
		    sprintf(filenames,"%s/sigma.txt",opt.path.c_str());
		}
		printf("%s\n",infilenames);
  		//FILE *outfile=fopen(filenames,"w");
//#if SUDAKOV>=1
		FILE *outfile=fopen(filenames,"w");
		//fclose(outfile);
		double val=0,x=0;
		double r;
		double arr[200]={0};
		for(int k=0; k<100;++k){
			x=X_MIN*pow(X_MAX/X_MIN,((double)k)/99);
			sigma.set_x(x);
			
//#pragma omp parallel private(r)
{
			//double k2,mu2=0;
//#pragma omp for schedule(dynamic)
			for(int j=0;j<100;++j){
				r=R_MIN*pow(R_MAX/R_MIN,((double)j)/99);
				printf("r= %.2e\n",r);
				arr[j]=sigma(x,r);
				
				printf("\033[1A\033[2K\r");
				if(arr[j]<0){
					printf("%.2e %.2e %.2e\n",x,r,arr[j]);
				}

			}
}
			for(int j=0;j<100;++j){
				r=R_MIN*pow(R_MAX/R_MIN,((double)j)/99);
				val=arr[j];
				fprintf(outfile ,"%.10e\t%.10e\t%.10e\n",log(x),log(r),val/sigpar[0] );
				
			}
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
