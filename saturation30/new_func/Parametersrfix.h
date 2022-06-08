/* GBW Starting parameter values, errors: sigma_0,lambda,x_0, C, mu2,g1 */
#if ((MODEL==0)||( MODEL==2 )) 
char * par_name[7]	= {"sigma_0",	"lambda",	"x_0",	"C", 	"r_max",	"g1",	"g2"};
double par_start[7]	= {	12.0,	0.29,	3.0e-4,	1.26,	0.5,	0.2,	0.8};
double par_error[7]	= {	1.0,	0.05,	0.1e-04,	0.01,	0.0,	0.01,	0.01 };

double   par_min[7]	= {	0.0,	0.00,	0.0,		0.1,	0.1 ,	0.0,	0.0};
double   par_max[7] 	= {	80.0,	1.00,	1.0,		10.0,	1.0,	2.0,	2.0};
#endif


