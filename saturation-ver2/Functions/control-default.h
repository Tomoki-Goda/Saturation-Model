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
	#define R_MIN 1.0e-5
#endif
#ifndef R_MAX
	#define R_MAX 30
#endif

#ifndef MU0
//if 1 mu02 is the fit parameter if 0 r_max is.
	#define MU0 1 
#endif

//////////////////IRREGULAR CONTROL//////////////////////
#ifndef  SATURATION
//if 0 use small r limit of model such that it does not have saturation 
	#define SATURATION  1
#endif



//////////////////////////////////////////////////////////////////
/////////////////////  system control ////////////////////////////
//////////////////////////////////////////////////////////////////
#ifndef PRINT_PROGRESS
	#define PRINT_PROGRESS 0
#endif

#ifndef N_SIMPS_R
//number of sampling for R integration, points are 2*N+1. //comparison with Fejer suggests it needs about 250 
	#define N_SIMPS_R 200 
#endif

#ifndef N_CHEB_R
//number of sampling for R integration, points are N divisible by 8
	#define N_CHEB_R 120 
#endif

#ifndef DGAUSS_PREC 
//precision of integration for adaptive gauss quadrature integration. or other methods
	#define DGAUSS_PREC 1.0e-4
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


#ifndef R_CHANGE_VAR
//use R=r/(1-r) for r integration. //seems to severely affect the quality...
	#define R_CHANGE_VAR 0
#endif

#ifndef NONLINEAR
//whether to use nonlinear transformation for fejer-curtis-clenshaw quadrature.
	#define NONLINEAR 0 
#endif

#ifndef STRATEGY
// Strategy for MIGRAD. read MINUIT documentation.
	#define STRATEGY 1
#endif

#ifndef FEJER
	#define FEJER 0
#endif

///////////////////   IRREGULAR CONTROL ///////////////////////
//// NOT TESTED, DISCONTINUED, ETC... CHECK WHEN CHANGED //////
///////////////////////////////////////////////////////////////
#ifndef NEW_DATA
	#define NEW_DATA 1 //1 mean only reading new hera
#endif

#ifndef Z_INTEGRATE
//integrate photon wave function over Z in advance
	#define Z_INTEGRATE 1
#endif

#ifndef TEST
// can be used for test... not properly defined...
	#define TEST 0
#endif

#ifndef SIMPS_GBS
//DISCONTINUED. use in-house simpson integration for GBS. 
	#define SIMPS_GBS 0
#endif

#ifndef SIMPS_Z_INT
//DISCONTINIED use in-house integration for z integration of photon psi.
	#define SIMPS_Z_INT 0
#endif

#ifndef THETA_OFF
//DISCONTINUED  . allow lower limit of sudakov to be larger than upper limit. 
	#define THETA_OFF 0
#endif

//////////////////////////////////////////////////////////////////
/////////// not to be chaged without a good reason ...///////////////////
//////////////////////////////////////////////////////////////////
//N_PAR is number of fit parameters
#if MODEL ==0 
	#define N_PAR 3
#elif MODEL==1
	#define N_PAR 5
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


//dipole sigma 
#if MODEL==0
#define SIGMA(r,x,Q2,par, sudpar)  sigma_gbw(r,x,Q2,par)
#elif MODEL==1
#define SIGMA(r,x,Q2,par, sudpar)  sigma_bgk(r,x,Q2,par)
#elif MODEL==2
#define SIGMA(r,x,Q2,par, sudpar)  sigma_gbs(r,x,Q2,par)  
#elif MODEL==22
#define SIGMA sigma_s 
#elif MODEL==3
#define SIGMA sigma_s
#endif

//dipole sigma to be combined with sudakov
#if SATURATION ==1
	#if (MODEL==3||MODEL==1)
		#define BASE_SIGMA sigma_bgk
	#elif (MODEL==22||MODEL==2||MODEL==0)
		#define BASE_SIGMA sigma_gbw
	#endif
#elif SATURATION==0
	#if (MODEL==3||MODEL==1)
		#define BASE_SIGMA sigma_bgk_ns
	#elif (MODEL==22||MODEL==2||MODEL==0)
		#define BASE_SIGMA sigma_gbw_ns
	#endif
#endif









