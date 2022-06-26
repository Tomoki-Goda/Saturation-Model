#include<math.h>
#include<stdio.h>
#include"control.h"
#include"../control-default.h"
#include"../constants.h"
//#include"./photon-wave-function-2.h"
#include<stdlib.h>


extern double psisq_z_int(double ,double ,int);

int main(int argc, char** argv){
	FILE *file;
	
	unsigned n=100;
	double val;
	double r;
	double Q2=atoi(argv[1]);
	char name[100];
	sprintf(name,"./photon-psi%s.txt",argv[1]);
	
	file=fopen(name,"w");
	for (unsigned i=0 ; i<=n;i++){
		r=pow(10,-2+2*((double)i)/n);
		val=pow(r,-2)*psisq_z_int(r,Q2,i);
		fprintf(file,"%f\t%f\n",r,val);
	}
	fclose(file);
}
