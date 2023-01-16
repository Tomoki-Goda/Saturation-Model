/* GBW Starting parameter values, errors: sigma_0,lambda,x_0, C, mu2,g1 */
#if ((MODEL==0)||( MODEL==2 )||( MODEL==22 )) 
std::string par_name[]	= {"sigma_0",	"lambda",	"x_0",		"mu102",	 "C2", 	"mu202",	"g1",	"g2"};
double par_start[]	= {3.0e+01,	3.3e-01,	2.0e+00,	4.0,		1.26,	2.0,		0.5,	0.5};
double par_error[]	= {	5.0,	0.05,		0.1,		1.0,		0.01,	0.1,		0.1,	0.1 };

double   par_min[]	= {	0.0,	0.00,		0.0,		1.0,		0.1,	1.0 ,		0.0,	0.0};
double   par_max[] 	= {	80.0,	1.00,		100.0,		10,		10.0,	10.0,		2.0,	2.0};
//#endif
#elif (MODEL==1)||(MODEL==3)
std::string  par_name[]	= {"sigma_0",	"A_g",	"lambda_g",	"mu102",	"C1", 	"mu02",		"C2",	"mu202",	"g1",	"g2"};
double par_start[]  	= {22.40,	1.0,	0.1,		4.0,		0.38, 	2.0,		1.26,	2.0,		0.1, 	0.1};  
double par_error[]  	= { 1.00,  	0.10,	0.05, 		1.0,		1.0,	0.1,		0.1,	0.1,		0.01,	0.01};
double   par_min[] 	= { 0.00,  	0.00,	-10.0,		1.0,		0.01, 	1.0,		0.5,	1.0,		0.0,	0.0};
double   par_max[] 	= {50.00,  	20.00,	10.0, 		10.0		10.0,	10.0,		10,	10.0,		2.0,	2.0};
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
int i=0,j=0,k=0;
#if (MODEL==0||MODEL==2||MODEL==22)
	sigpar[j++]=2.56819*(PREC)(par[i++]);//mb to GeV
	sigpar[j++]=(PREC)(par[i++]);
	sigpar[j++]=(PREC)(par[i++])*1.0e-4;

#else
	sigpar[j++]=2.56819*(PREC)(par[i++]);
	sigpar[j++]=(PREC)(par[i++]);
	//sigpar[1]=pow(fabs((PREC)par[1])*1.0e-4,par[2]);
	sigpar[j++]=(PREC)(par[i++]);
#endif

#if MU02==0
	sigpar[j++]=par[i++];
#else 
	++i;
#endif
	//printf("SIGMA: %.2e %.2e %.2e ",sigpar[0],sigpar[1],sigpar[2]);
#if(MODEL==1||MODEL==3)
	sigpar[j++]=(PREC)(par[i++]);
	sigpar[j++]=(PREC)(par[i++]);//sqrt(fabs(sigpar[3]/par[4]));//rmax^2= C/mu02
	//printf(" %.2e %.2e ",sigpar[3],sigpar[4]);
#endif
////////////////////////////SUDPAR////////////////////////////////
#if (MODEL==22||MODEL==2)
	sudpar[k++]=(PREC)(par[i++]);
	sudpar[k++]=(PREC)(par[i++]);
	//printf("\tSUDAKOV: %.2e %.2e ",sudpar[0],sudpar[1]);	
	#if (SUDAKOV==2)
		sudpar[k++]=(PREC)(par[i++]);
		sudpar[k++]=(PREC)(par[i++]);
		//printf("%.2e %.2e ",sudpar[2],sudpar[3]);
	#endif

///////////////////////////////////////////////////////
#elif (MODEL==3)
	#if INDEPENDENT_C==1
		sudpar[k++]=(PREC)(par[i++]) ;
	#else 
		sudpar[k++]=(PREC)(par[(i++)-2]);
	#endif

	#if INDEPENDENT_RMAX==1
		sudpar[k++]=(PREC)(par[i++]);
	#else
		sudpar[k++]=(PREC)(par[(i++)-2]);//mu02 is shared
	#endif
	//printf("%.2e %.2e ",sudpar[0],sudpar[1]);

//////////////////////////////////////////////////////
#if (SUDAKOV==2)
	sudpar[k++]=(PREC)(par[i++]);
	sudpar[k++]=(PREC)(par[i++]);
	//printf("%.2e %.2e ",sudpar[2],sudpar[3]);
#endif
#endif
	//printf("\n");
	if(i!=N_PAR){
		printf("parameter number mismatch counted %d, N_PAR= %d\n",i,N_PAR);
		exit(0);
	}
	return 0;
}
