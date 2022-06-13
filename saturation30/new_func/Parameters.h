/* GBW Starting parameter values, errors: sigma_0,lambda,x_0, C, mu2,g1 */
#if ((MODEL==0)||( MODEL==2 )) 
char * par_name[7]	= {"sigma_0",	"lambda",	"x_0",	"C", 	"r_max",	"g1",	"g2"};
double par_start[7]	= {	12.0,	0.29,	3.0e-4,	1.26,	0.5,	0.2,	0.8};
double par_error[7]	= {	1.0,	0.05,	0.1e-04,	0.01,	0.01,	0.01,	0.01 };

double   par_min[7]	= {	0.0,	0.00,	0.0,		0.1,	0.1 ,	0.0,	0.0};
double   par_max[7] 	= {	80.0,	1.00,	1.0,		10.0,	1.0,	2.0,	2.0};

#elif (MODEL==1)
char * par_name[5]	= {"sigma_0",	"A_g",	"lambda_g",	"C", 	"r_max"};
double par_start[5]  	= {22.40,	1.35,	0.079,		0.38, 	1.73}; 
double par_error[5]  	= { 1.00,  	0.10,	0.05, 		0.05, 	0.40};
double   par_min[5] 	= { 0.00,  	0.00,	0.01, 		0.01, 	0.50};
double   par_max[5] 	= {30.00,  	5.00,	0.20, 		1.00, 	5.00};
#endif


