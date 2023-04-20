//////////////////////////////////////////////////////////////////
/////////////////////// regular control //////////////////////////
//////////////////////////////////////////////////////////////////
//////////////     MODEL AND FIXED PARAMETERS             ////////
//////////////////////////////////////////////////////////////////
#ifndef MODEL 
	#define MODEL 0
#endif

#ifndef SUDAKOV 
	#define SUDAKOV 0
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
#ifndef X_DATA_MAX 
	#define X_DATA_MAX 1.0e-2
#endif
#ifndef Q2_MAX 
	#define Q2_MAX 1.0e+1
#endif

////// upper and lower cut off of r ////////
//Keep in mind r and k limits can be related.  
//generally R_MAX should be larger than 1/sqrt(KT2_MIN) etc
////////////////////////////////////////
#ifndef R_MIN 
	#define R_MIN 1.0e-6
#endif
#ifndef R_MAX
	#define R_MAX 1e+5
#endif


#ifndef X_MIN 
	#define X_MIN 1.0e-8
#endif
#ifndef X_MAX 
	#define X_MAX 1.0e+00
#endif
#ifndef KT2_MIN 
	#define KT2_MIN 1.0e-6
#endif

//////////////////////////////////////////////////////////////////
/////////////////////  system control ////////////////////////////
//////////////////////////////////////////////////////////////////

#ifndef WW 
	#define WW 0
#endif 
#ifndef ADJOINT 
	#if WW==1
		#define ADJOINT 1
	#else 
		#define ADJOINT 0
	#endif
#endif

#ifndef PRINT_PROGRESS
	#define PRINT_PROGRESS 0
#endif

#ifndef N_CHEB_R
//number of sampling for R integration, points are N divisible by 8
	#define N_CHEB_R 100 
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


#ifndef MU02
	#define MU02 1
#endif

#ifndef GLUON_APPROX
	#define GLUON_APPROX 1
#endif


#ifndef FREEZE_QS2 
	#define FREEZE_QS2 0
#endif

#ifndef ADD_END
	#define ADD_END 0
#endif


#ifndef R_CHANGE_VAR
//use R=r/(1-r) for r integration. 
	#define R_CHANGE_VAR 0
#endif

#ifndef USE_RESULT//use result.txt. Select 2 to round up the result.
	#define USE_RESULT 0
#endif

#ifndef R_FORMULA//dipole factorization
	#define R_FORMULA 0
#endif
#ifndef SIGMA_APPROX
	#define SIGMA_APPROX -1
#endif

#ifndef THRESHOLD//threshold factor (1-x)^a
	#define THRESHOLD 0
#endif
#ifndef NS//non-pert. sudakov like factor in the fourier-hankel trans-integral. see dipole-gluon.hh
	#define NS 0
#endif
#ifndef VARIANT //variants of models. see r-formula.h
	#define VARIANT 0
#endif

#ifndef N_CHEB//variants of models. see r-formula.h
	#define N_CHEB 25
#endif
#ifndef CHEB_D//variants of models. see r-formula.h
	#define CHEB_D 1
#endif

#ifndef SECTOR_MAX
	#define SECTOR_MAX 200
#endif
//////////////////////////////////////////////////////////////////
/////////// not to be chaged without a good reason ...///////////////////
//////////////////////////////////////////////////////////////////
//N_PAR is number of fit parameters
#if MODEL ==0 
	#define N_PAR 5
#elif MODEL==1
	#define N_PAR 7
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







