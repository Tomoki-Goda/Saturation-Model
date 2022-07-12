#include<stdio.h>
#include<math.h>
#include<time.h>
#include<cfortran.h>
#include<minuit.h>

#include"control.h"
#include"control-default.h"
#include"constants.h"

#include"./Parameters.h"

#define MAXN 600
extern double SIGMA_PREC;
extern unsigned N_DATA;
extern int N_SIMPS;

extern int load_data(void);
extern void generate_psi_set(void);
//extern void fcn(int , double , double, double ,unsigned ,void (*)(void) );
extern void fcn(int* npar, double grad[], double*fcnval, double *par,unsigned* iflag,void (*dum)(void) );

double arglist[10];

void log_printf(FILE* file,char* line){
	if(file!=stdout){
		fprintf(file,"%s",line);
	}
#if (PRINT_PROGRESS==1)
	fprintf(stdout,"%s",line);
#endif
}


FILE* log_file;
FILE * out_file;
//////////////////////////////////////////////////////////////////////////
//////////////////////////  main  ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]){
#if (R_FIX==1)
	par_error[4]=0.0; 
#endif
	int n_data;
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
	
	
	/////////////////   for results   //////////////////////
	char name[11];//apparently thats the max length minuit accepts
	double res;
	//double res_par[N_PAR];
	double error;
	double dum3;
	int dum4;

	int error_flag = 0;
	///////////////////// Initial /////////////////////////
	
	//load_data( Q2_data, x_data, y_data, w_data, cs_data, err_data,  N_DATA);
	n_data=load_data();
	//printf("%d\n",N_DATA);
	printf("%d\n",n_data);

	
	printf("----------------------------------  Starting  --------------------------------\n");
	printf("Model ID:  %d  \t Q2_up: %.1e \t x_up: %.1e \t  Sudakov: %d\n", MODEL, Q2_MAX,X_MAX, SUDAKOV);
	printf("R_FIX: %d \t N_PAR %d                                     \n",R_FIX,N_PAR );
	printf("L %.2e S %.2e C %.2e B %.2e\n",MASS_L2,MASS_S2,MASS_C2,MASS_B2 );
	printf(" STAR %d\n", STAR );
	printf("Gauss eps: %.2e\t Simps N: %d \t \n", DGAUSS_PREC,N_SIMPS_R);
	printf("-----------------------------------------------------------------------------\n");
	
	
//////////////////////////     Initialize Minuit     ////////////////////////////////////////
	MNINIT(5,6,7);
	/* Parameters definition */
	for(unsigned i=0;i<N_PAR;i++){
		MNPARM(i+1,par_name[i], par_start[i],par_error[i],par_min[i],par_max[i],error_flag);
	}
	
	
	char command[100];
////////////////////////  Minimalization procedure MIGRAD //////////////////////////////////
	//sprintf(command , "SET ERRORDEF 3");
	//sprintf(command , "MIGRAD 500 %d %.3e",10*(N_PAR*N_PAR),100.0);
	//sprintf(command , "SET EPSMACHINE 1.0e-5");
	//MNCOMD(fcn,command,error_flag,0);



	N_SIMPS=(int)(N_SIMPS_R*2.0/4.0);
	SIGMA_PREC=DGAUSS_PREC*50;
	generate_psi_set();
	
	/* Set strategy to STRategy from main.h, 0-fast, 1-default, 2-precise */
	MNCOMD(fcn,"SET STR 0", error_flag,0);
	
	sprintf(command , "SIMPLEX %d 2.5D-1 ",10*(N_PAR*N_PAR));//get first digit right
	MNCOMD(fcn,command,error_flag,0);
	
	sprintf(command , "SET LIMITS");//Having limits seriously deteriorates performance
	MNCOMD(fcn,command,error_flag,0);
	
	//SIGMA_PREC=DGAUSS_PREC*10;
	
	sprintf(command , "MIGRAD %d %.3e",100*(N_PAR*N_PAR),5.0);
	MNCOMD(fcn,command,error_flag,0);

#if (PRINT_PROGRESS==1)	
	sprintf(command , "SET PRINTOUT 3");
	MNCOMD(fcn,command,error_flag,0);
#endif
////////////////////////////////   FINAL RUN    ///////////////////////////////////////	
	SIGMA_PREC=DGAUSS_PREC;
	N_SIMPS=N_SIMPS_R;
	generate_psi_set();
	//sprintf(command,"SET ERRORDEF %.3e", 0.25); 
	//MNCOMD(fcn,command,error_flag,0);
	
	sprintf(command,"SET STR %d", STRATEGY); 
	MNCOMD(fcn,command,error_flag,0);
	
	
	for(int i=0; i<3;i++){
	
		sprintf(command , "MIGRAD %d %.3e",200*(N_PAR*N_PAR),2.5);
		MNCOMD(fcn,command,error_flag,0);
		if(error_flag==0){
			break;
		}
	}
	//sprintf(command , "IMPROVE %d ",50*(N_PAR*N_PAR) );
	//MNCOMD(fcn,command,error_flag,0);
	//if(error_flag!=0){
	//	MNCOMD(fcn,command,error_flag,0);
	//}
	
	time_measure-=clock();
/////////////////////////////////////SAVE RESULTS////////////////////////////////
	
	char outline[500];
	sprintf(outline,"Q_up\t%.0f\n",Q2_MAX);
	log_printf(out_file,outline);
	double res_par[N_PAR];
	for(unsigned i=0;i<N_PAR;i++){
		MNPOUT(i+1,name,res,error,dum3,dum3,dum4);
		sprintf(outline,"%s\t%.4e\t%.4e\n",name,res,error);
		log_printf(out_file,outline);
		*(res_par+i)=res;
	}
	
	MNSTAT(res,error,dum3,dum4,dum4,dum4);
	
	sprintf(outline,"chisq\t%.4e\t%.4e\n",res,error);
	log_printf(out_file,outline);
	sprintf(outline,"n_data-n_par\t%d\n",n_data-N_PAR);
	log_printf(out_file,outline);
	sprintf(outline,"chisq/dof\t%.3e\n",res);///(n_data-N_PAR));
	log_printf(out_file,outline);
	
	
	sprintf(outline,"\n\n");
	log_printf(out_file,outline);
	sprintf(outline,"In %.2e minutes\n", -((double)time_measure)/(60*CLOCKS_PER_SEC) );
	log_printf(out_file,outline);
	
	sprintf(outline,"Error Flag %d \n", error_flag );
	log_printf(out_file,outline);
	
	
	fclose(out_file);
	fclose(log_file);
	
	printf("*************************  End   **************************\n");
	printf("Model ID:  %d  \t Q2_up: %.1e \t x_up: %.1e \t  Sudakov: %d\n", MODEL, Q2_MAX,X_MAX, SUDAKOV);
	printf("R_FIX: %d \t                                               \n",R_FIX );
	printf("L %.2e S %.2e C %.2e B %.2e\n",MASS_L2,MASS_S2,MASS_C2,MASS_B2 );
	printf(" STAR %d\n", STAR );
	printf("Gauss eps: %.2e\t Simps N: %d \t \n", DGAUSS_PREC,N_SIMPS_R);
	printf("chisq/dof\t%.3e\n",res);
	printf("****************************************************************\n");
	return(0);
 
}


