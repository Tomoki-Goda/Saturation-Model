
//////////////////////////////////////////////////////////////////
/////////////////////  system control ////////////////////////////
//////////////////////////////////////////////////////////////////
#ifndef Z_INTEGRATE
	#define Z_INTEGRATE 1
#endif

#ifndef TEST	
	#define TEST 0
#endif

#ifndef SIMPS_GBS 
	#define SIMPS_GBS 0
#endif

#ifndef N_SIMPS_R
	#define N_SIMPS_R 250
#endif

#ifndef R_CHANGE_VAR
	#define R_CHANGE_VAR 1
#endif

#ifndef NEW_DATA
	#define NEW_DATA 1 //1 mean only reading new hera
#endif

#ifndef PRINT_PROGRESS
	#define PRINT_PROGRESS 0
#endif

#ifndef DGAUSS_PREC
	#define DGAUSS_PREC 0.001
#endif
//////////////////////////////////////////////////////////////////
/////////// not to be chaged without a good reason ...///////////////////
//////////////////////////////////////////////////////////////////

#if MODEL ==0 
	#define N_PAR 3
#elif MODEL==1
	#define N_PAR 5
#elif MODEL==2
	#if SUDAKOV==0
		#define N_PAR 3
	#elif SUDAKOV==1
		#define N_PAR 5
	#elif SUDAKOV==2
		#define N_PAR 7
	#endif
#endif



#if MODEL==0
#define SIGMA sigma_gbw
#elif MODEL==1
#define SIGMA sigma_bgk
#elif MODEL==2
#define SIGMA sigma_gbs
#endif

#define STRATEGY 0.0
