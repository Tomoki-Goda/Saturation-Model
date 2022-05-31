//////////////////////////////////////////////////////////////////
/////////////////////// regular control //////////////////////////
//////////////////////////////////////////////////////////////////
#define X_MAX 1.0e-2
#define Q2_MAX 5.0e+1

#define MODEL 2
#define FLAVOUR 2
#define SUDAKOV 2

/* GBW Starting parameter values, errors: sigma_0,lambda,x_0, C, mu2,g1 */
#if ((MODEL==0)||( MODEL==2 )) 
char * par_name[7]	= {"sigma_0",	"lambda",	"x_0",	"C", 	"r_max",	"g1",	"g2"};
double par_start[7]	= {	12.0,	0.29,	3.0e-4,	1.26,	0.5,	0.2,	0.8};
double par_error[7]	= {	1.0,	0.05,	0.1e-04,	0.01,	0.0,	0.01,	0.01 };

double   par_min[7]	= {	0.0,	0.00,	0.0,		0.01,	0.0 ,	0.0,	0.0};
double   par_max[7] 	= {	80.0,	1.00,	1.0,		20.0,	2.0,	2.0,	2.0};
#endif






#define MASS_L2 0.0
#define MASS_S2 0.0
#define MASS_C2 1.96
#define MASS_B2 21.16


//////////////////////////////////////////////////////////////////
/////////////////////  system control ////////////////////////////
//////////////////////////////////////////////////////////////////
#define Z_INTEGRATE 1 
#define TEST 0
#define SIMPS_GBS 1

#define N_SIMPS_R 200

//////////////////////////////////////////////////////////////////
/////////// not to be chaged without a good reason ...///////////////////
//////////////////////////////////////////////////////////////////

#if MODEL ==0 
#define N_PAR 3
#elif MODEL==1
#define N_PAR 5
#elif MODEL==2
#define N_PAR 7
#endif


#if MODEL==0
#define SIGMA sigma_gbw
#elif MODEL==1
#define SIGMA sigma_bgk
#elif MODEL==2
#define SIGMA sigma_gbs
#endif

#define STRATEGY 0.0
