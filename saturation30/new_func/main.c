#include<stdio.h>
#include<math.h>
#include<time.h>

#include"./control_tmp.h"
#include"./constants.h"

#include"simpson-integral.h"
#include"dipole-cross-section.h"
#include"photon-wave-function-2.h"

#include"cfortran.h"
#include"../minuit.h"
#include"./read-and-fit.h"

//#define TEST 2
#define MAXN 600

static double x_data[MAXN]={0};
static double y_data[MAXN]={0};
static double w_data[MAXN]={0};
static double Q2_data[MAXN]={0};
static double cs_data[MAXN]={0};
static double err_data[MAXN]={0};
static unsigned N_DATA;

double arglist[10];





//////////////////////////////////////////////////////////////////////////
//////////////////////////  main  ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]){
	clock_t time_measure=clock();
	
	FILE * out_file;
	if(argc>1){
		//FILE * out_file= fopen("./results.txt","w");
		out_file= fopen(argv[1],"w");
	}else{
		//out_file=stdout;
		out_file= fopen("./result.txt","w");
	}

	/////////////////results//////////////////////
	char name[11];//apparently thats the max length minuit accepts
	double res;
	double error;
	double dum3;
	int dum4;

	int error_flag = 0;
	
	load_data( Q2_data, x_data, y_data, w_data, cs_data, err_data,  N_DATA);
	generate_psi_set();
	
	
#if TEST==1
	double resval=0.0;	
	double grad[7];
	fcn(7, grad, &resval , par_start,0, &dum_func );
	printf("fcn returns %f,\n",resval); 
	return(0);
#endif
	
	/* Initialize Minuit */
	MNINIT(5,6,7);
	/* Parameters definition */
	for(unsigned i=0;i<N_PAR;i++){
		MNPARM(i+1,par_name[i], par_start[i],par_error[i],par_min[i],par_max[i],error_flag);
	}
	arglist[0] = STRATEGY;
	/* Set strategy to STRategy from main.h, 0-fast, 1-default, 2-precise */
	MNEXCM(fcn,"SET STR",arglist, 1,error_flag,0);
	
	/* Minimalization procedure MIGRAD */
	MNEXCM(fcn,"MIGRAD",0,0,error_flag,0);
	
	////////////////////SAVE RESULTS////////////////////////////////
	
	for(unsigned i=0;i<N_PAR;i++){
	
		MNPOUT(i+1,name,res,error,dum3,dum3,dum4);
		fprintf(out_file,"%s\t%.4e\t%.4e\n",name,res,error);
	}
	
	
	
	
	MNSTAT(res,error,dum3,dum4,dum4,dum4);
	fprintf(out_file,"chisq\t%.4e\t%.4e\n",res,error);
	fprintf(out_file,"N_DATA\t%d\t%.3e\n",N_DATA,res/(N_DATA-N_PAR));
	
	
	time_measure-=clock();
	
	fprintf(out_file,"\n\n");
	fprintf(out_file,"In %.2e minutes\n", -((double)time_measure)/(60*CLOCKS_PER_SEC) );
	fclose(out_file);
	return(0);
 
}


