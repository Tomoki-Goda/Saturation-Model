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

#include"minuit-run.h"

#define TEST_M 0
//////////////////////////////////////////////////////////////////////////
//////////////////////////  main  ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]){
	char outline[100];
	int error_flag;
	double corr[10];//10 is just enough to hold all par;
	double err_mat[10*10];
	double dum;
	
#if (R_FIX==1)
	par_error[4]=0.0; 
#endif
//	int n_data;
	clock_t time_measure=clock();
	//FILE* logfile;
	//FILE * out_file;
	char resultfile[500];
	char logfile[500];
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
	
	/////////////Initialize etc//////////////////
	MINUIT_INIT();
#if TEST_M==1
	CHECK_COV(&fcn);
	return 0;
#endif
	RUN_MINUIT(&fcn);
	//int error_flag;
	//MNCOMD(*fcn,"SET LIMITS",error_flag,0);
	//MNCOMD(*fcn,"MINIMIZE 1000 1.0D00 " ,error_flag,0);
	MNCOMD(*fcn,"CALLFCN 3",error_flag,0);	
	time_measure-=clock();
	
	
	sprintf(outline,"\n\nCorrelation:\n");
	log_printf(log_file,outline);
	for(int i=0;i<N_PAR;i++){
		corr[i]=0;
		MNERRS(i+1,dum,dum,dum,corr[i]);
		sprintf(outline,"%.4f  ",corr[i]);
		log_printf(log_file,outline);
				
	}
	sprintf(outline,"\n\n");
	log_printf(log_file,outline);
	
	MNEMAT(*err_mat,N_PAR);
	for(int i=0;i<N_PAR;i++){
		for(int j=0;j<N_PAR;j++){
			sprintf(outline,"%.4e\t",err_mat[i*N_PAR+j]);
			log_printf(log_file,outline);
		}
		sprintf(outline,"\n\n");
		log_printf(log_file,outline);
	}
/////////////////////////////////////SAVE RESULTS////////////////////////////////
	out_file=fopen(resultfile,"w");
	SAVE_RESULT(out_file);
	
	sprintf(outline,"\n\n");
	log_printf(out_file,outline);
	sprintf(outline,"In %.2e minutes\n", -((double)time_measure)/(60*CLOCKS_PER_SEC) );
	log_printf(out_file,outline);
	
	fclose(out_file);
	fclose(log_file);
	return 0;	
}


