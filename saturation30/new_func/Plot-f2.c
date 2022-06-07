#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include"./control_tmp.h"
#include"./control-default.h"
#include"./constants.h"

#include"dipole-cross-section.h"
#include"./photon-wave-function-2.h"
#include"./simpson-integral.h"

double f2_integrand(double r, double ** par){
	double xm;
	double x = *(*par);
	double Q2= *(*(par)+1);
	double *param;
	param=(*(par+1));
	double value=0.0;
	for(unsigned fl=0;fl< NF-1;fl++){
		//printf("%d",fl);
		xm=mod_x(x,Q2,fl);
		value+=psisq_z_int(r, Q2, fl)* SIGMA(r,xm,Q2, param)/r ;
		
	}
	//printf("x: %f\t Q2: %f\t%f\t%f\t%f\tresult :%f \n", x,Q2,*(param),*(param+1),*(param+2),value);
	return(value);
}

void generate_data_set(double *par,double xpow ,char* datafile){
	//csarray is counterpart of CS_DATA ...
	//double integral[N_DATA];
	double gevtofm = 0.1973;
	double point;
	double r,xm ;
	FILE* file=fopen(datafile,"w");
	printf("generate_data_set\n");
	unsigned point_n=50;
	unsigned point_n2=2;//take two points 
	double *param[2];//just big enough
	double var[10];
	(*param)=var;
	
	double grad;
	
	(*(param+1))=par;
	
	if(file==NULL){
		printf("generate_data_set::file error\n");
	}
	double res[2], err[2];
	//printf("s_0: %f\t l: %f\t x_0: %f \n", *(*(param+1)),(*(*(param+1)+1) ),(*(*(param+1)+2) ), res);
	
	for(unsigned i=0; i<point_n;i++){
		
		Q2=pow(10,-2+5*((double)i)/point_n);
		*(var+1) =Q2;//Q2
		for(unsigned j=0;j<2;j++){
			*(res+j)=0;
			(*(var))= pow( 10,xpow+0.1*j );
			
			simpson1dA(&f2_integrand,param,1.0e-6,30,100,res+j,err+j);
			//printf("NF:%f\n",NF);	
		}
		grad=(log((*res)) - log( (*(res+1))) )/ (0.1 *log(10) );
		
		fprintf(file,"%lf\t%lf\t%lf\n",Q2,*res, grad);	
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
		if( strcmp(argv[2*i-1],"-x")==0 ){
			x=log(strtod(argv[2*i],&end)) /log(10);
		}else if( strcmp(argv[2*i-1],"-out")==0 ){
			sprintf(resultfilename,"%s",argv[2*i]);
		}else if( strcmp(argv[2*i-1],"-in" ) ==0){
			sprintf(parfilename,"%s" ,argv[2*i]);
		}else{
			printf("Please Use flag, -x, -out , -in\n\n");
		}
			
	}
	printf("x=%f\n",x);
	
	
	printf("%s\n%s\n",parfilename,resultfilename);
	FILE* parfile=fopen(parfilename,"r");
	
	if(parfile==NULL){
		printf("Plot-F2.c file error\n");
		return 1;
	}
	
/////////////////////read result data //////////////////////////////
	fscanf(parfile,"%s\t%f",name,&dum);//line for Qup
	
	for(unsigned i=0;i<N_PAR;i++){
		fscanf(parfile,"%s\t%lf\t%f\n",name,param+i,&dum);
		fprintf(stdout,"%s \t\t %.3e  \n",name ,*(param+i));
	}
	fclose(parfile);
	generate_data_set(param, x ,resultfilename);
	
	return 0;

}


