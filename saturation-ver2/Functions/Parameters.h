/* GBW Starting parameter values, errors: sigma_0,lambda,x_0, C, mu2,g1 */
#if ((MODEL==0)||( MODEL==2 )||( MODEL==22 )) 
char * par_name[8]	= {"sigma_0",	"lambda",	"x_0",	"C", 	"r_max",	"g1",	"g2"};
double par_start[8]	= {	23.0,	0.29,	3.0e-4,	1.26,	0.4,	0.5,	0.5};
double par_error[8]	= {	1.0,	0.05,	0.1e-04,	0.01,	0.01,	0.1,	0.1 };

double   par_min[8]	= {	0.0,	0.00,	0.0,		0.1,	0.1 ,	0.0,	0.0};
double   par_max[8] 	= {	80.0,	1.00,	1.0,		10.0,	1.0,	2.0,	2.0};

#elif (MODEL==1)
char * par_name[8]	= {"sigma_0",	"A_g",	"lambda_g",	"C", 	"r_max"	};
double par_start[8]  	= {22.40,	1.35,	0.079,		0.38, 	0.4		}; 
double par_error[8]  	= { 1.00,  	0.10,	0.05, 		0.1, 	0.01		};
double   par_min[8] 	= { 0.00,  	0.00,	0.01, 		0.1, 	0.1		};
double   par_max[8] 	= {50.00,  	5.00,	0.20, 		10.00,	1.00		};

#elif (MODEL==3)
#if (INDEPENDENT_C==0)
char * par_name[8]	= {"sigma_0",	"A_g",	"lambda_g",	"C", 	"r_max",	"C2",	"g1",	"g2"};
double par_start[8]  	= {22.40,	1.35,	0.079,		0.38, 	0.4,		1.0,	0.1, 	0.1}; 
double par_error[8]  	= { 1.00,  	0.10,	0.05, 		0.01, 	0.01,		0.0,	0.01,	0.01};
double   par_min[8] 	= { 0.00,  	0.00,	0.01, 		0.01, 	0.01,		1.0,	0.0,	0.0};
double   par_max[8] 	= {50.00,  	5.00,	0.20, 		10.00,	1.00,		50,	2.0,	2.0};
#else 
char * par_name[8]	= {"sigma_0",	"A_g",	"lambda_g",	"C", 	"r_max",	"C2",	"g1",	"g2"};
double par_start[8]  	= {22.40,	1.0,	0.05,		0.5, 	0.5,		1.0,	0.1, 	0.1};  
//double par_start[8]  	= {22.30,	1.64,	0.0101,	0.20, 	0.30,		20.0,	0.5, 	0.5}; 
double par_error[8]  	= { 1.00,  	0.10,	0.01, 		0.01, 	0.01,		0.1,	0.01,	0.01};
double   par_min[8] 	= { 0.00,  	0.00,	0.00, 		0.01, 	0.01,		0.5,	0.0,	0.0};
double   par_max[8] 	= {50.00,  	5.00,	0.20, 		10.00,	1.00,		10,	2.0,	2.0};
#endif
#endif


