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

int main(){
	FILE* outfile1=fopen("./dp-plots.m","w");
	FILE* outfile2=fopen("./f2-plots.m","w");

	double par[7]={1.0,0.24,0.5e-4,1.26,4.0,0.2,0.8  };

	double Qs2=10;
	double x=1.0e-3;
	
	double r, rmin, rmax;
	rmax=10.0;
	rmin=0.001;
	unsigned point_n=1000;
	double pos, val1, val2;
	double step=( (double)(rmax-rmin) )/point_n;
	fprintf(outfile1,"{");
	fprintf(outfile2,"{");
	
	clock_t total;
	clock_t t;
	double Q2;
	for(unsigned i=0;i<=point_n;i++){
		t=clock();
		pos=pow(10.0 , -3+((double)6*i)/1000);
		
		Q2=pow(10.0 , -3+((double)6*i)/1000);
		
		val1=sigma_gbw(pos, mod_x(x,Qs2,'l'), Qs2, par );
		
		val2=sigma_DIS(x,Q2,1.0,par);
		
		
		t-=clock();
		total-=t;

		fprintf(outfile1, "{%f,%f}",pos, val1);
		(i==point_n)?(fprintf(outfile1,"}")):fprintf(outfile1,",") ;
		
		fprintf(outfile2, "{%f,%f}",Q2, val2);
		(i==point_n)?(fprintf(outfile2,"}")):fprintf(outfile2,",") ;
	}
	printf("%f seconds\n",((float)total )/CLOCKS_PER_SEC);
	
	fclose(outfile1);
	fclose(outfile2);
	
	system("math -i <\"plotdp.wls\"");
	system("math -i <\"plotf2.wls\"");

}
