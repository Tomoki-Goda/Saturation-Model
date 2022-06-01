#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include"./control_tmp.h"
#include"./control-default.h"
#include"./constants.h"

#include"dipole-cross-section.h"

void generate_data_set(double *par, double Q2,double x ,char* datafile){
	//csarray is counterpart of CS_DATA ...
	//double integral[N_DATA];
	double gevtofm = 0.1973;
	double point;
	double r,xm;
	FILE* file=fopen(datafile,"w");
	printf("generate_data_set\n");
	unsigned point_n=10000;
	
	for(unsigned i=0; i<point_n;i++){
		point=0;
		r =pow(10,-2+3*((double)i)/point_n);
		//printf("%f\n",r);
		for(unsigned fl=0;fl<(NF-1);fl++){
				xm=mod_x(x, Q2,fl );
				point+=  ( SIGMA(r,xm,Q2, par) );
							
		
		}
		fprintf(file,"%f\t%f\n",gevtofm*r,point/((NF-1)*(*par)));		
	}
	fclose(file);
}


int main(int argc, char ** argv){
	char parfilename[100];
	char resultfilename[100];
	char name[20];
	float par;
	double param[10];//just any number large enough
	float dum;
	//char dumc[100];
	
	if(argc>1){
		//FILE * out_file= fopen("./results.txt","w");
		sprintf(parfilename,"%s/result.txt",argv[1]);
		sprintf(resultfilename,"%s/plot.txt",argv[1]);
	}else{
		sprintf(parfilename,"./result.txt");
		sprintf(resultfilename,"./plot.txt");		
	}
	FILE* parfile=fopen(parfilename,"r");
	
	fscanf(parfile,"%s\t%f\n",name,&par);//first line is Qup
	printf("%s\t%.0f\n",name,par);

	for(unsigned i=0;i<N_PAR;i++){
		//printf("%d\n",i);
		fscanf(parfile,"%s\t%lf\t%f\n",name,param+i,&dum);
		fprintf(stdout,"%s \t\t %.3e  \n",name ,*(param+i));
	}
	fclose(parfile);
	generate_data_set(param, 100.0,0.001 ,resultfilename);
	
	return 0;

}


