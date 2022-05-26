#include<stdio.h>
#include<math.h>
#include<time.h>

#include"./control.h"
#include"./constants.h"

#include"simpson-integral.h"
#include"dipole-cross-section.h"
#include"photon-wave-function.h"

#include"DIS-cross-section.h"

#include"cfortran.h"
#include"../minuit.h"
#include"./read-and-fit.h"
//#include<file.h>

//FILE* outfile=fopen("./plots.m","w");
#define PLOT 2

int main(){
	double par[7]={23.0,0.28,3.0e-4,1.26,4.0,0.2,0.8  };

	
	double test_z,test_r,test_Q2, test_val, test_x;
	test_z=0.0001;
	test_Q2=1000.0;
	test_x=1.0e-4;
	unsigned point_n=50;
#if PLOT==0
/////////////////////////////////////////////////////////////////
	FILE* testfile=fopen("./dp.txt","w");
	for(unsigned i =0;i<point_n;i++){
		test_r =pow(10,-2+3*((double)i)/point_n);
		test_val=SIGMA(test_r,mod_x(test_x,test_Q2,'l'), test_Q2, par) ;
		test_val*=1/(par[0]);
		printf("%.2e %.3e\n",test_r,test_val);
		fprintf(testfile,"%f\t%f\n",test_r,test_val);
	}
	fclose(testfile);
#elif PLOT==1
///////////////////////////////////////////////////////////////
	FILE* testfile1=fopen("./pwf.txt","w");
	
	for(unsigned i =0;i<point_n;i++){
		test_r =pow(10,-3+4*((double)i)/point_n);
		test_val=psisq_f ( test_r,test_z,  test_Q2, 'l'); 
		fprintf(testfile1,"%f\t%f\n",test_r,test_val);
		
	}
	fclose(testfile1);
#else
/////////////////////////////////////////////////////////////////
	FILE* outfile2=fopen("./f2-plots.txt","w");
	
	double r, rmin, rmax;
	rmax=10.0;
	rmin=0.001;
	
	double pos, val1, val2;
	double step=( (double)(rmax-rmin) )/point_n;
	//fprintf(outfile1,"{");
	//fprintf(outfile2,"{");
	
	clock_t total;
	clock_t t;
	
	for(unsigned i=0;i<=point_n;i++){
		t=clock();		
		test_Q2=pow(10.0 , -3+((double)6*i)/point_n);
				
		val2= sigma_DIS(test_x,test_Q2,1.0,par);
		t-=clock();
		total-=t;

		fprintf(outfile2, "%f\t%f",test_Q2, val2);
		(i==point_n)?(fprintf(outfile2," ")):fprintf(outfile2,"\n") ;
	}
	printf("%f seconds\n",((float)total )/CLOCKS_PER_SEC);
	
	//fclose(outfile1);
	fclose(outfile2);
#endif
	
}
