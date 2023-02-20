//////////////////////////////////////////////////////////////////
/////////////////////// regular control //////////////////////////
//////////////////////////////////////////////////////////////////
//////////////     MODEL AND FIXED PARAMETERS             ////////
//////////////////////////////////////////////////////////////////
#ifndef MODEL 
	#define MODEL 2
#endif
#ifndef FLAVOUR 
	#define FLAVOUR 2
#endif
#ifndef SUDAKOV 
	#define SUDAKOV 2
#endif

#ifndef MASS_L2
	#define MASS_L2 0.0196
#endif
#ifndef MASS_S2 
	#define MASS_S2 0.0196
#endif
#ifndef MASS_C2 
//	#define MASS_C2 1.96
	#define MASS_C2 1.69
#endif
#ifndef MASS_B2
	#define MASS_B2 21.16
#endif
//////////////SELECTION OF DATA////////////
#ifndef X_MAX 
	#define X_MAX 1.0e-2
#endif
#ifndef Q2_MAX 
	#define Q2_MAX 1.0e+1
#endif

////// upper and lower cut off of r ////////
#ifndef R_MIN
	#define R_MIN 1.0e-6
#endif
#ifndef R_MAX
	#define R_MAX 1e+5
#endif

//////////////////////////////////////////////////////////////////
/////////////////////  system control ////////////////////////////
//////////////////////////////////////////////////////////////////
#ifndef PRINT_PROGRESS
	#define PRINT_PROGRESS 0
#endif

#ifndef N_CHEB_R
//number of sampling for R integration, points are N divisible by 8
	#define N_CHEB_R 250 
#endif

#ifndef STAR
//star presctiption for r. 0 is Collins=Soper type, 1 is Golec-Biernat=Sapeta type
	#define STAR 1
#endif

#ifndef R_FIX
//treat r_max or mu02 as constant.
	#define R_FIX 0
#endif

#ifndef INDEPENDENT_C
//treat two C in model3 C/r independently
	#define INDEPENDENT_C 1
#endif

#ifndef INDEPENDENT_RMAX
//treat r_max or mu02 in model 3 independently
	#define INDEPENDENT_RMAX 0
#endif



#ifndef LAPLACIAN
	#define LAPLACIAN 0
#endif
#ifndef IBP
	#define IBP 0
#endif

#ifndef ALPHA_RUN
	#define ALPHA_RUN 0 
#endif
#ifndef MODX
	#define MODX 0
#endif
#ifndef PHI
	#define PHI 0
#endif

#ifndef SCATTER
	#define SCATTER 0
#endif

#ifndef MU02
	#define MU02 1
#endif

#ifndef GLUON_APPROX
	#define GLUON_APPROX 1
#endif

#ifndef HANKEL
	#define HANKEL 0
#endif
#ifndef FREEZE_QS2 
	#define FREEZE_QS2 0
#endif

#ifndef ADD_END
	#define ADD_END 0
#endif
#ifndef GBW_APPROX
	#define GBW_APPROX 0
#endif

#ifndef R_CHANGE_VAR
//use R=r/(1-r) for r integration. 
	#define R_CHANGE_VAR 0
#endif

#ifndef USE_RESULT
	#define USE_RESULT 0
#endif


//////////////////////////////////////////////////////////////////
/////////// not to be chaged without a good reason ...///////////////////
//////////////////////////////////////////////////////////////////
//N_PAR is number of fit parameters
#if MODEL ==0 
	#define N_PAR 4
#elif MODEL==1
	#define N_PAR 6
#elif (MODEL==2||MODEL==22)
	#if SUDAKOV==0
		#define N_PAR 3
	#elif SUDAKOV==1
		#define N_PAR 5
	#elif SUDAKOV==2
		#define N_PAR 7
	#endif
#elif MODEL==3
	#if SUDAKOV==0
		#define N_PAR 5
	#elif SUDAKOV==1
		#if ((INDEPENDENT_C==1)&&(INDEPENDENT_RMAX==1))
			#define N_PAR 7
		#elif ( ((INDEPENDENT_C==0)&&(INDEPENDENT_RMAX==1))||( (INDEPENDENT_C==1)&&(INDEPENDENT_RMAX==0)) )
			#define N_PAR 6
		#else
			#define N_PAR 5
		#endif

	#elif SUDAKOV==2
		#if ((INDEPENDENT_C==1)&&(INDEPENDENT_RMAX==1))
			#define N_PAR 9
		#elif ( ((INDEPENDENT_C==0)&&(INDEPENDENT_RMAX==1))||( (INDEPENDENT_C==1)&&(INDEPENDENT_RMAX==0)) )
			#define N_PAR 8
		#else
			#define N_PAR 7
		#endif
	#endif
#endif







