#define MAXN 600

extern double SIGMA_PREC;
static unsigned N_DATA;
extern int N_SIMPS;
extern int N_CHEB;
extern int N_OFF;//number of arameters fixed

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

int MINUIT_INIT(){ 
	int error_flag = 0;
	///////////////////// Initial /////////////////////////
	N_DATA=load_data(); 
	printf("%d\n",N_DATA);
	
	printf("----------------------------------  Starting  --------------------------------\n");
	printf("Model ID:  %d  \t Q2_up: %.1e \t x_up: %.1e \t  Sudakov: %d\n", MODEL, Q2_MAX,X_MAX, SUDAKOV);
	printf("R_FIX: %d \t N_PAR %d                                     \n",R_FIX,N_PAR );
	printf("L %.2e S %.2e C %.2e B %.2e\n",MASS_L2,MASS_S2,MASS_C2,MASS_B2 );
	printf("STAR %d R_CHANGE_VAR %d\n", STAR,R_CHANGE_VAR );
	printf("Gauss eps: %.2e\t Simps N: %d CHEB N: %d \t \n", DGAUSS_PREC,N_SIMPS_R, N_CHEB_R);
	printf("-----------------------------------------------------------------------------\n");
	
	
//////////////////////////     Initialize Minuit     ////////////////////////////////////////
	MNINIT(5,6,7);
	/* Parameters definition */
	int skip=0;
	for(unsigned i=0;(i-skip)<N_PAR;i++){
#if MODEL==3
///////////////
//parameters  may be shared bet. BGK & Sud, see dipole-cross-section.h parameters() for how they are organized.
///////////////
#if INDEPENDENT_C==0
		if(i==5){
			skip++;
			continue;
		}
#endif
#if INDEPENDENT_RMAX==0
		if(i==6){
			skip++;
			continue;
		}
#endif
#endif

		MNPARM(i+1-skip,par_name[i], par_start[i],par_error[i],par_min[i],par_max[i],error_flag);
	}
	
	return error_flag;	
}



/*
int RUN_MINUIT(void(*fcn)(int* , double*, double*, double *,unsigned*,void (*)(void) ) ){
	char command[100];
	int error_flag, istat, nvpar, npar;
	double val, edm, up;

	N_SIMPS=(int)(N_SIMPS_R*3.0/5.0);
	SIGMA_PREC=DGAUSS_PREC*10;
	generate_psi_set();
	
	MNCOMD(*fcn, "SET PRINTOUT 3",error_flag,0);
#if (MODEL==3||MODEL==1)		
	MNCOMD(*fcn, "FIX 3",error_flag,0);
#if SUDAKOV>=1			
	MNCOMD(*fcn,"FIX 5",error_flag,0);
#endif

#elif (MODEL==2||MODEL==22)
	MNCOMD(*fcn, "FIX 3",error_flag,0);
#if SUDAKOV>=1
	MNCOMD(*fcn, "FIX 5",error_flag,0);
#endif
#endif

	MNCOMD(*fcn, "SIMPLEX 1000 2",error_flag,0);
/////////////////////////////////////////////////////////////////////
	
	N_SIMPS=(int)(N_SIMPS_R);
	SIGMA_PREC=DGAUSS_PREC;
	generate_psi_set();
	
	MNCOMD(*fcn,"SET LIMITS",error_flag,0);

	MNCOMD(*fcn,"HESSE",error_flag,0);
	MNCOMD(*fcn, "MIGRAD 1000 1",error_flag,0);

	MNSTAT(val,edm, up, nvpar,npar,istat);
	printf("\n***********  FIRST RUN  **********\n");
	printf("ISTAT= %d     FCN/DOF=%.3e     EDM=%.3e\n",istat, val/(N_DATA-N_PAR),edm);	
	printf("************************************\n");
	
	//N_SIMPS=N_SIMPS_R;
	//SIGMA_PREC=DGAUSS_PREC;
	//generate_psi_set();
	
#if (MODEL==3||MODEL==1)			
	MNCOMD(*fcn, "RELEASE 3",error_flag,0);
#if SUDAKOV>=1
	MNCOMD(*fcn, "RELEASE 5",error_flag,0);
#endif	
#elif (MODEL==2||MODEL==22)			
	MNCOMD(*fcn, "RELEASE 3",error_flag,0);
#if SUDAKOV>=1
	MNCOMD(*fcn, "RELEASE 5",error_flag,0);
#endif
#endif
	MNCOMD(*fcn,"HESSE",error_flag,0);
	MNCOMD(*fcn, "MIGRAD 1000 1",error_flag, 0);

	MNSTAT(val,edm, up, nvpar,npar,istat);
	printf("\n*********   FINAL RUN   ***********\n");
	printf("ISTAT= %d     FCN/DOF=%.3e     EDM=%.3e\n",istat, val/(N_DATA-N_PAR),edm);	
	printf("************************************\n");

	return istat;
}

*/
int SAVE_RESULT(FILE* outfile){
	char outline[500];
	
	char name[11];
	double res_par[N_PAR];
	double res, error,dum3;
	int dum4, istat;	

	sprintf(outline,"Q_up\t%.0f\n",Q2_MAX);
	log_printf(out_file,outline);
	
	for(unsigned i=0;i<N_PAR;i++){
		MNPOUT(i+1,name,res,error,dum3,dum3,dum4);
		sprintf(outline,"%s\t%.4e\t%.4e\n",name,res,error);
		log_printf(out_file,outline);
		*(res_par+i)=res;
	}
	
	

	
	MNSTAT(res,error,dum3,dum4,dum4,istat);
	sprintf(outline,"chisq\t%.4e\t%.4e\n",res,error);
	log_printf(out_file,outline);
	sprintf(outline,"n_data\t%d\n",N_DATA);
	log_printf(out_file,outline);
	sprintf(outline,"chisq/dof\t%.3e\n",res/(N_DATA-(N_PAR-N_OFF)));
	log_printf(out_file,outline);
	
	
	
	sprintf(outline,"Error Flag %d \n", istat );
	log_printf(out_file,outline);

	printf("*************************  End   **************************\n");
	printf("Model ID:  %d  \t Q2_up: %.1e \t x_up: %.1e \t  Sudakov: %d\n", MODEL, Q2_MAX,X_MAX, SUDAKOV);
	printf("R_FIX: %d \t                                               \n",R_FIX );
	printf("L %.2e S %.2e C %.2e B %.2e\n",MASS_L2,MASS_S2,MASS_C2,MASS_B2 );
	printf("STAR %d R_CHANGE_VAR %d\n", STAR,R_CHANGE_VAR );
	printf("Gauss eps: %.2e\t Simps N: %d CHEB N: %d \t \n", DGAUSS_PREC,N_SIMPS_R, N_CHEB_R);
	printf("chisq/dof\t%.3e \nError Matrix %d\n",res/(N_DATA-N_PAR)  , istat);
	printf("****************************************************************\n");
	return(0);
 
}

