#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//#include<time.h>
#include<string.h>

//#include"control.h"
//#include"../control-default.h"
//#include"../constants.h"

int read_options(int argc, char ** argv, double* param,double* x,double* Q2,char* outfile ){
	char parfilename[100]="";
	char name[20]="";
	float par;
	//double param[10];//just any number large enough
	float dum;
	char* end;
	
	//for(int i=0;i<argc;i++){
	//	printf("%s\t ",argv[i]);
	//}
	//printf("\n");

	for(unsigned i=1;i<=argc/2;i++){
		if( strcmp(argv[2*i-1],"-x")==0 ){
		//	*x=log(strtod(argv[2*i],&end)) /log(10);
			*x=strtod(argv[2*i],&end) ;
			*x=pow(10.0,-(*x));
		}else if( strcmp(argv[2*i-1],"-k")==0 ){//use x for k //this should be fixed in the future
			*x=strtof(argv[2*i],&end) ;
		}else if( strcmp(argv[2*i-1],"-Q2" ) ==0){
			*Q2=strtof(argv[2*i],&end);
		}else if( strcmp(argv[2*i-1],"-W" ) ==0){
			*Q2=strtof(argv[2*i],&end);
		}else if( strcmp(argv[2*i-1],"-out")==0 ){
			sprintf(outfile,"%s",argv[2*i]);
		}else if( strcmp(argv[2*i-1],"-in" ) ==0){
			sprintf(parfilename,"%s" ,argv[2*i]);
		}else if( strcmp(argv[2*i-1],"-fl" ) ==0){
			*Q2=((float)atoi(argv[2*i]))+ 0.5;
			printf("plot fl= %f",*Q2);
		}else{
			printf("Please Use flag, -x, -out , -in\n\n");
		}
			
	}
	//printf("x=%.2e\n",*x);
	//printf("Q2=%.2e\n",*Q2);
		
	printf("%s\n%s\n",parfilename,outfile);
	FILE* parfile=fopen(parfilename,"r");
	
	if(parfile==NULL){
		printf("plot.c file error\n");
		return 1;
	}
	
/////////////////////read result data //////////////////////////////
	fscanf(parfile,"%s\t%f",name,&dum);//line for Qup
	
	//for(unsigned i=0;i<N_PAR;i++){
	for(int i=0;i<10/*any large enough*/;i++){
		fscanf(parfile,"%s\t%lf\t%f\n",name,param+i,&dum);
		if(strcmp(name,"chisq")==0){
			break;
		}
		fprintf(stdout,"%s \t\t %.3e  \n",name ,*(param+i));
	}
	fclose(parfile);
	
	

	//generate_data_set(param, x ,resultfilename);
	
	return 0;

}

int plot(double (*func)(double ,double**), double* var ,int varlen,  double ** par ,FILE* file ){
	//printf("plot %d\n",varlen);
	//FILE* file=fopen(file_name,"w");
//	if(file==NULL){
//		printf(" plot:: file error. %s.\n", file_name);
//	}	
	double val;
	
	for(int i=0 ; i<varlen; i++){
		val=(*func)(*(var+i),par);	
		//fprintf(file,"%f\t%f\n",*(var+i),val);
		fprintf(file,"%.5e\t%.5e\n",*(var+i),val);
		//printf("%.2e\t%.2e\n", *(var+i),val);

	}
	
	//fclose(file);
	return 0;
}
