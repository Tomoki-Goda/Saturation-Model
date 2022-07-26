/* GBW Starting parameter values, errors: sigma_0,lambda,x_0, C, mu2,g1 */
#if ((MODEL==0)||( MODEL==2 )||( MODEL==22 )) 
#if MU0==1
char * par_name[8]	= {"sigma_0",	"lambda",	"x_0",		"C", 	"mu02",		"g1",	"g2"};
double par_start[8]	= {	23.0,	0.29,		3.0,		1.26,	5.0,		0.5,	0.5};
double par_error[8]	= {	1.0,	0.05,		0.01,		0.01,	0.1,		0.1,	0.1 };

double   par_min[8]	= {	0.0,	0.00,		0.0,		0.1,	1.0 ,		0.0,	0.0};
double   par_max[8] 	= {	80.0,	1.00,		100.0,		10.0,	10.0,		2.0,	2.0};
#else
char * par_name[8]	= {"sigma_0",	"lambda",	"x_0",		"C", 	"r_max",	"g1",	"g2"};
double par_start[8]	= {	23.0,	0.29,		3.0,		1.26,	0.5,		0.5,	0.5};
double par_error[8]	= {	1.0,	0.05,		0.01,		0.01,	0.01,		0.1,	0.1 };

double   par_min[8]	= {	0.0,	0.00,		0.0,		0.1,	0.1 ,		0.0,	0.0};
double   par_max[8] 	= {	80.0,	1.00,		100.0,		10.0,	1.0,		2.0,	2.0};
#endif

#elif (MODEL==1)
#if MU0==1
char * par_name[8]	= {"sigma_0",	"A_g",	"lambda_g",	"C", 	"mu02"	};
double par_start[8]  	= {22.40,	1.35,	0.1,		0.38, 	2.0		}; 
double par_error[8]  	= { 1.00,  	0.10,	0.05, 		0.1, 	0.1		};
double   par_min[8] 	= { 0.00,  	0.00,	-10.0,		0.1, 	1.0		};
double   par_max[8] 	= {50.00,  	5.00,	10.0, 		10.00,	10.00		};
#else
char * par_name[8]	= {"sigma_0",	"A_g",	"lambda_g",	"C", 	"r_max"	};
double par_start[8]  	= {22.40,	1.35,	0.1,		0.38, 	0.4		}; 
double par_error[8]  	= { 1.00,  	0.10,	0.05, 		0.1, 	0.01		};
double   par_min[8] 	= { 0.00,  	0.00,	-10.0,		0.1, 	0.1		};
double   par_max[8] 	= {50.00,  	5.00,	10.0, 		10.00,	1.00		};
#endif

#elif (MODEL==3)
#if MU0==1
char * par_name[]	= {"sigma_0",	"A_g",	"lambda_g",	"C", 	"mu02",	"C2",	"mu202",	"g1",	"g2"};
double par_start[]  	= {22.40,	1.0,	0.1,		0.38, 	2.0,		1.26,	5.0,		0.1, 	0.1};  
double par_error[]  	= { 1.00,  	0.10,	0.05, 		0.01, 	0.1,		0.1,	0.1,		0.01,	0.01};
double   par_min[] 	= { 0.00,  	0.00,	-10.0,		0.01, 	1.0,		0.5,	1.0,		0.0,	0.0};
double   par_max[] 	= {50.00,  	5.00,	10.0, 		10.00,	10.0,		10,	10.0,		2.0,	2.0};
#else
char * par_name[]	= {"sigma_0",	"A_g",	"lambda_g",	"C", 	"r_max",	"C2",	"r_max2",	"g1",	"g2"};
double par_start[]  	= {22.40,	1.0,	0.1,		0.38, 	0.4,		1.26,	0.5,		0.1, 	0.1};  
double par_error[]  	= { 1.00,  	0.10,	0.05, 		0.01, 	0.01,		0.1,	0.01,		0.01,	0.01};
double   par_min[] 	= { 0.00,  	0.00,	-10.0,		0.01, 	0.01,		0.5,	0.01,		0.0,	0.0};
double   par_max[] 	= {50.00,  	5.00,	10.0, 		10.00,	1.00,		10,	1.00,		2.0,	2.0};
#endif
#endif


