//////////////////////////////////////////////////////////////////
/////////////////////// regular control //////////////////////////
//////////////////////////////////////////////////////////////////
#define X_MAX 1.0e-2
#define Q2_MAX 1.0e+2

#define MODEL 2
#define FLAVOUR 2
#define SUDAKOV 2

/* GBW Starting parameter values, errors: sigma_0,lambda,x_0, C, mu2,g1 */
#if ((MODEL==0)||( MODEL==2 )) 
char * par_name[7]	= {"sigma_0",	"lambda",	"x_0",	"C", 	"r_max",	"g1",	"g2"};
double par_start[7]	= {	24.0,	0.29,	3.0e-4,	1.26,	0.5,	0.2,	0.8};
double par_error[7]	= {	1.0,	0.05,	0.1e-04,	0.01,	0.0,	0.01,	0.01 };

double   par_min[7]	= {	0.0,	0.00,	0.0,		0.01,	0.0 ,	0.0,	0.0};
double   par_max[7] 	= {	80.0,	1.00,	1.0,		20.0,	2.0,	2.0,	2.0};
#endif






#define MASS_L2 0.0196
#define MASS_S2 0.0196
#define MASS_C2 1.96
#define MASS_B2 21.16


