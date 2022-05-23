#include<stdio.h>
#include<math.h>
#include<time.h>

#include"./constants.h"
#define MODEL 2 
#define FLAVOUR 2


#include"simpson-integral.h"
#include"dipole-cross-section.h"
#include"photon-wave-function.h"

#include"DIS-cross-section.h"




int main(){
	fprintf(stdout,"started\n");

	FILE* file=fopen("./sigma.m", "w");
	fprintf(file, "{");
	double x=1.0e-3;
		
	double par[7]={1.0,0.24,0.5e-4,1.26,4.0,0.2,0.8  };

	double res;
	double q2;
	clock_t time;
	for(unsigned i=0;i<100;i++){
		q2=1.0e-3+10*i;
		time=clock();
		res=sigma_DIS(x,q2,1,par);
		time-=clock();
		fprintf(file, "{%f,%f},",q2,res);
		fprintf(stdout, "{%.3e,%.3e} : %.3e ,\n",q2,res, -((double)time)/CLOCKS_PER_SEC);

	}
	fprintf(file,"\b}");
	fclose(file);


	return 0;
}

