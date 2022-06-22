#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include"../control_tmp.h"
#include"../control-default.h"
#include"../constants.h"

#include"../gluon-chebyshev.h"
#include"../dipole-cross-section.h"

#if (MODEL==3||MODEL==22||MODEL==2)
#include"../sudakov.h"
#endif


//Syntax is ./main -in <input dir > -out <output dir > -Q2 <Q2 > -x <x >

void generate_data_set(double *par, double Q2,double x ,char* datafile){
	//csarray is counterpart of CS_DATA ...
	//double integral[N_DATA];
	double gevtofm = 0.1973;
	double point;
	double r,xm;
	FILE* file=fopen(datafile,"w");
	printf("generate_data_set\n");
	unsigned point_n=100;
	if(file==NULL){
		printf("generate_data_set::file error\n");
		
	}
	
	
	for(int i=0; i<point_n;i++){
		point=0;
		r =pow(10,-2+3*((double)i)/point_n);
		//printf("%f\n",r);
		//printf("x: %f\nQ2: %f\n", x,Q2);
		for(unsigned fl=0;fl<(NF-1);fl++){
		//for(int fl=0;fl<1;fl++){
				xm=mod_x(x, Q2,fl );
				//xm=x;
				point+=  ( SIGMA(r,xm,Q2, par) );
							
		
		}
		fprintf(file,"%f\t%f\n",gevtofm*r,point/( (NF-1) * (*par)));		
	}
	fclose(file);
}


int main(int argc, char ** argv){
	char parfilename[100]="";
	char resultfilename[100]="";
	char name[20]="";
	float par;
	double param[10];//just any number large enough
	float dum;
	//char dumc[100];
	double Q2=100;
	double x=0.001;
	char* end;
	
	for(unsigned i=1;i<=argc/2;i++){
		if( strcmp(argv[2*i-1],"-Q2")==0 ){
			Q2=strtod(argv[2*i],&end);
		}else if( strcmp(argv[2*i-1],"-x")==0 ){
			x=strtod(argv[2*i],&end);
		}else if( strcmp(argv[2*i-1],"-out")==0 ){
			sprintf(resultfilename,"%s",argv[2*i]);
		}else if( strcmp(argv[2*i-1],"-in" ) ==0){
			sprintf(parfilename,"%s" ,argv[2*i]);
		}else{
			printf("Please Use flag, -Q2, -x, -out , -in\n\n");
		}
			
	}
	printf("Q2= %f, x=%f\n",Q2,x);
	
	
	printf("%s\n%s\n",parfilename,resultfilename);
	FILE* parfile=fopen(parfilename,"r");
	
	if(parfile==NULL){
		printf("Plot-DP.c file error\n");
		return 1;
	}
	
	fscanf(parfile,"%s\t%f",name,&dum);//line for Qup
	
	for(unsigned i=0;i<N_PAR;i++){
		fscanf(parfile,"%s\t%lf\t%f\n",name,param+i,&dum);
		fprintf(stdout,"%s \t\t %.3e  \n",name ,*(param+i));
	}
	fclose(parfile);
	
#if (MODEL==1||MODEL==3)	
	approx_xg(param+1);//generate chebyshev coefficients
#endif
	generate_data_set(param, Q2,x ,resultfilename);
	
	return 0;

}


