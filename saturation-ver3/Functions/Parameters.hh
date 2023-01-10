/* GBW Starting parameter values, errors: sigma_0,lambda,x_0, C, mu2,g1 */
#if ((MODEL==0)||( MODEL==2 )||( MODEL==22 )) 
#if MU0==1
std::string par_name[8]	= {"sigma_0",	"lambda",	"x_0",		"C2", 	"mu202",	"g1",	"g2"};
double par_start[8]	= {3.0e+01,	3.0e-01,	2.0e+00,	1.26,	2.0,		0.5,	0.5};
double par_error[8]	= {	1.0,	0.05,		0.01,		0.01,	0.1,		0.1,	0.1 };

double   par_min[8]	= {	0.0,	0.00,		0.0,		0.1,	1.0 ,		0.0,	0.0};
double   par_max[8] 	= {	80.0,	1.00,		100.0,		10.0,	10.0,		2.0,	2.0};
#else
char * par_name[8]	= {"sigma_0",	"lambda",	"x_0",		"C2", 	"r_max",	"g1",	"g2"};
double par_start[8]	= {	23.0,	0.29,		3.0,		1.26,	0.5,		0.5,	0.5};
double par_error[8]	= {	1.0,	0.05,		0.01,		0.01,	0.01,		0.1,	0.1 };

double   par_min[8]	= {	0.0,	0.00,		0.0,		0.1,	0.1 ,		0.0,	0.0};
double   par_max[8] 	= {	80.0,	1.00,		100.0,		10.0,	1.0,		2.0,	2.0};
#endif

#elif (MODEL==1)
#if MU0==1
std::string par_name[8]	= {"sigma_0",	"A_g",	"lambda_g",	"C", 	"mu02"	};
double par_start[8]  	= {22.40,	1.35,	0.1,		0.38, 	2.0		}; 
double par_error[8]  	= { 1.00,  	0.10,	0.05, 		0.1, 	0.1		};
double   par_min[8] 	= { 0.00,  	0.00,	-10.0,		0.1, 	1.0		};
double   par_max[8] 	= {50.00,  	20.00,	10.0, 		10.00,	10.00		};
#else
std::string  par_name[8]	= {"sigma_0",	"A_g",	"lambda_g",	"C", 	"r_max"	};
double par_start[8]  	= {22.40,	1.35,	0.1,		0.38, 	0.4		}; 
double par_error[8]  	= { 1.00,  	0.10,	0.05, 		0.1, 	0.01		};
double   par_min[8] 	= { 0.00,  	0.00,	-10.0,		0.1, 	0.1		};
double   par_max[8] 	= {50.00,  	20.00,	10.0, 		10.00,	1.00		};
#endif

#elif (MODEL==3)
#if MU0==1
std::string  par_name[]	= {"sigma_0",	"A_g",	"lambda_g",	"C", 	"mu02",	"C2",	"mu202",	"g1",	"g2"};
double par_start[]  	= {22.40,	1.0,	0.1,		0.38, 	2.0,		1.26,	2.0,		0.1, 	0.1};  
double par_error[]  	= { 1.00,  	0.10,	0.05, 		0.01, 	0.1,		0.1,	0.1,		0.01,	0.01};
double   par_min[] 	= { 0.00,  	0.00,	-10.0,		0.01, 	1.0,		0.5,	1.0,		0.0,	0.0};
double   par_max[] 	= {50.00,  	20.00,	10.0, 		10.00,	10.0,		10,	10.0,		2.0,	2.0};
#else
std::string  par_name[]	= {"sigma_0",	"A_g",	"lambda_g",	"C", 	"r_max",	"C2",	"r_max2",	"g1",	"g2"};
double par_start[]  	= {22.40,	1.0,	0.1,		0.38, 	0.4,		1.26,	0.5,		0.1, 	0.1};  
double par_error[]  	= { 1.00,  	0.10,	0.05, 		0.01, 	0.01,		0.1,	0.01,		0.01,	0.01};
double   par_min[] 	= { 0.00,  	0.00,	-10.0,		0.01, 	0.01,		0.5,	0.01,		0.0,	0.0};
double   par_max[] 	= {50.00,  	20.00,	10.0, 		10.00,	1.00,		10,	1.00,		2.0,	2.0};
#endif
#endif

int parameter(const std::vector<double>& par,PREC(& sigpar)[],PREC(& sudpar)[]){
///////////////////
//Sigpar are as we all know it, parameters for dipole sigma.
//sudpar are {C , r_max, g1, g2} but parameters may be given in terms of mu02 (as in BGK), and C and r_max may be that of BGK.
//That's why this is so messy...
////////////////////
	//for(int i=0;i<N_PAR;i++){
	//	printf("%.2e ",par[i]);
	//}
	//printf("\n");
#if (MODEL==0||MODEL==2||MODEL==22)
	sigpar[0]=2.56819*(PREC)(par[0]);//mb to GeV
	sigpar[1]=(PREC)(par[1]);
	sigpar[2]=(PREC)(par[2])*1.0e-4;
#else
	sigpar[0]=2.56819*(PREC)(par[0]);
	sigpar[1]=(PREC)(par[1]);
	//sigpar[1]=pow(fabs((PREC)par[1])*1.0e-4,par[2]);
	sigpar[2]=(PREC)(par[2]);
#endif

	
	
	//printf("SIGMA: %.2e %.2e %.2e ",sigpar[0],sigpar[1],sigpar[2]);
#if(MODEL==1||MODEL==3)
	sigpar[3]=(PREC)(par[3]);
	#if MU0==0
		sigpar[4]=sigpar[3]/((PREC)(par[4])*(PREC)(par[4]));
	#else
		sigpar[4]=(PREC)(par[4]);//sqrt(fabs(sigpar[3]/par[4]));//rmax^2= C/mu02
	#endif
	//printf(" %.2e %.2e ",sigpar[3],sigpar[4]);
#endif
////////////////////////////SUDPAR////////////////////////////////
#if (MODEL==22||MODEL==2)
		sudpar[0]=(PREC)(par[3]);
	#if MU0==0
		sudpar[1]=sudpar[0]/((PREC)(par[4])*(PREC)(par[4]));
	#else
		sudpar[1]=(PREC)(par[4]);
	#endif
	//printf("\tSUDAKOV: %.2e %.2e ",sudpar[0],sudpar[1]);	
#if (SUDAKOV==2)
	sudpar[2]=(PREC)(par[5]);
	sudpar[3]=(PREC)(par[6]);
	//printf("%.2e %.2e ",sudpar[2],sudpar[3]);
#endif

///////////////////////////////////////////////////////
#elif (MODEL==3)
	#if INDEPENDENT_C==1
		sudpar[0]=(PREC)(par[5]) ;
	#else 
		sudpar[0]=(PREC)(par[3]);
	#endif

	#if MU0==0 //if rmax is fit parameter
		#if INDEPENDENT_RMAX==1
			sudpar[1]=sudpar[0]/((PREC)(par[6])*(PREC)(par[6]));
		#else
			sudpar[1]=sudpar[0]/((PREC)(par[4])*(PREC)(par[4]));
		#endif

	#else //if mu02 is the fit parameter
		#if INDEPENDENT_RMAX==1
			sudpar[1]=(PREC)(par[6]);
		#else
			sudpar[1]=(PREC)(par[4]);//mu02 is shared
		#endif
	#endif
	//printf("%.2e %.2e ",sudpar[0],sudpar[1]);

//////////////////////////////////////////////////////
#if (SUDAKOV==2)
	sudpar[2]=(PREC)(par[7]);
	sudpar[3]=(PREC)(par[8]);
	//printf("%.2e %.2e ",sudpar[2],sudpar[3]);
#endif
#endif
	//printf("\n");
	return 0;
}
