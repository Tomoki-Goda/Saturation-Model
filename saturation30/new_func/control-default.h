//////////////////////////////////////////////////////////////////
/////////////////////// regular control //////////////////////////
//////////////////////////////////////////////////////////////////
#ifndef X_MAX 
	#define X_MAX 1.0e-2
#endif
#ifndef Q2_MAX 
	#define Q2_MAX 1.0e+1
#endif

#ifndef MODEL 
	#define MODEL 2
#endif
#ifndef FLAVOUR 
	#define FLAVOUR 2
#endif
#ifndef SUDAKOV 
	#define SUDAKOV 2
#endif

/* GBW Starting parameter values, errors: sigma_0,lambda,x_0, C, mu2,g1 */
#ifndef MASS_L2
	#define MASS_L2 0.0196
#endif
#ifndef MASS_S2 
	#define MASS_S2 0.0196
#endif
#ifndef MASS_C2 
	#define MASS_C2 1.96
#endif
#ifndef MASS_B2
	#define MASS_B2 21.16
#endif

#ifndef R_FIX
	#define R_FIX 0
#endif
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
	#define N_SIMPS_R 150
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

#ifndef SIMPS_Z_INT
	#define SIMPS_Z_INT 0
#endif

#ifndef THETA_OFF
	#define THETA_OFF 0
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
