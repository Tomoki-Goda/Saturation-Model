#include<stdio.h>
#include<math.h>
#include<time.h>
#include<cfortran.h>
#include<minuit.h>

#include"control.h"
#include"control-default.h"
#include"constants.h"

#include"./Parameters.h"
#include"./main.h"

#define TEST_M 1
//////////////////////////////////////////////////////////////////////////
//////////////////////////  main  ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]){
	char outline[100];
#if (R_FIX==1)
	par_error[4]=0.0; 
#endif
//	int n_data;
	clock_t time_measure=clock();
	//FILE* logfile;
	//FILE * out_file;
	char resultfile[100];
	char logfile[100];
	//char datafile[100];
	if(argc>1){
		//FILE * out_file= fopen("./results.txt","w");
		sprintf(resultfile,"%s/result.txt",argv[1]);
		sprintf(logfile,"%s/log.txt",argv[1]);
		//sprintf(resultfile,"%s/plot.txt",argv[1]);
	}else{
		strcpy(resultfile,"./result.txt");
		strcpy(logfile,"./log.txt");
		//strcpy(resultfile,"./plot.txt");	
	}
	log_file=fopen(logfile,"w");
	out_file=fopen(resultfile,"w");
	
	/////////////Initialize etc//////////////////
	MINUIT_INIT();
#if TEST_M==1
	CHECK_COV(&fcn);
	return 0;
#endif
	RUN_MINUIT(&fcn);
	
	time_measure-=clock();
/////////////////////////////////////SAVE RESULTS////////////////////////////////
	SAVE_RESULT(out_file);
	
	sprintf(outline,"\n\n");
	log_printf(out_file,outline);
	sprintf(outline,"In %.2e minutes\n", -((double)time_measure)/(60*CLOCKS_PER_SEC) );
	log_printf(out_file,outline);
	
	fclose(out_file);
	fclose(log_file);
	return 0;	
}


