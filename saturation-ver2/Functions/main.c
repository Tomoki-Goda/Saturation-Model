#include<stdio.h>
#include<math.h>
#include<time.h>

//#include"./control_tmp.h"
#include"control.h"
#include"control-default.h"
#include"constants.h"



//#include"gluon-chebyshev.h"
//#include"simpson-integral.h"
//#include"dipole-cross-section.h"

//#if (MODEL==3||MODEL==22||MODEL==2)
//#include"sudakov.c"
//#endif


//#include"photon-wave-function-2.h"


//in the directory /usr/include ??
#include"cfortran.h"
#include"minuit.h"


#include"./Parameters.h"




//#else 
//#include"./Parametersrfix.h"
//#endif

//#if MODEL==1
//#include"../gluons.h"
//#include"../chebyshev.h"
//#/include"../chebyshev3.h"
//#endif

//#define TEST 2
#define MAXN 600
extern int load_data(void);
extern void generate_psi_set(void);
//extern void fcn(int , double , double, double ,unsigned ,void (*)(void) );
extern void fcn(int npar, double grad[], double*fcnval, double *par,unsigned iflag,void (*dum)(void) );

//static double x_data[MAXN]={0};
//static double y_data[MAXN]={0};
//static double w_data[MAXN]={0};
//static double Q2_data[MAXN]={0};
//static double cs_data[MAXN]={0};
//static double err_data[MAXN]={0};
//extern unsigned N_DATA;

double arglist[10];

void log_printf(FILE* file,char* line){
	if(file!=stdout){
		fprintf(file,"%s",line);
	}
#if PRINT_PROGRESS==1
	fprintf(stdout,"%s",line);
#endif
}


FILE* log_file;
FILE * out_file;
//////////////////////////////////////////////////////////////////////////
//////////////////////////  main  ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]){
#if R_FIX==1
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
	generate_psi_set();
	
	printf("*************************  Starting  **************************\n");
	printf("Model ID:  %d  \t Q2_up: %.1e \t x_up: %.1e \t  Sudakov: %d\n", MODEL, Q2_MAX,X_MAX, SUDAKOV);
	printf("R_FIX: %d \t N_PAR %d                                     \n",R_FIX,N_PAR );
	printf("L %.2e S %.2e C %.2e B %.2e\n",MASS_L2,MASS_S2,MASS_C2,MASS_B2 );
	printf(" STAR %d\n", STAR );
	printf("Gauss eps: %.2e\t Simps N: %d \t \n", DGAUSS_PREC,N_SIMPS_R);
	printf("****************************************************************\n");
	
	
	
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
	
	time_measure-=clock();
	////////////////////SAVE RESULTS////////////////////////////////
	
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
	sprintf(outline,"chisq/dof\t%.3e\n",res/(n_data-N_PAR));
	log_printf(out_file,outline);
	
	
	sprintf(outline,"\n\n");
	log_printf(out_file,outline);
	sprintf(outline,"In %.2e minutes\n", -((double)time_measure)/(60*CLOCKS_PER_SEC) );
	log_printf(out_file,outline);
	fclose(out_file);
	fclose(log_file);
	
	printf("*************************  Koniec   **************************\n");
	printf("Model ID:  %d  \t Q2_up: %.1e \t x_up: %.1e \t  Sudakov: %d\n", MODEL, Q2_MAX,X_MAX, SUDAKOV);
	printf("R_FIX: %d \t                                               \n",R_FIX );
	printf("L %.2e S %.2e C %.2e B %.2e\n",MASS_L2,MASS_S2,MASS_C2,MASS_B2 );
	printf(" STAR %d\n", STAR );
	printf("Gauss eps: %.2e\t Simps N: %d \t \n", DGAUSS_PREC,N_SIMPS_R);
	printf("chisq/dof\t%.3e\n",res/(n_data-N_PAR));
	printf("****************************************************************\n");
	return(0);
 
}


