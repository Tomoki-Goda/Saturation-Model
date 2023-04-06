/* GBW Starting parameter values, errors: sigma_0,lambda,x_0, C, mu2,g1 */
#if ((MODEL==0)||( MODEL==2 )||( MODEL==22 )) 
std::string par_name[]	= {"sigma_0",	"lambda",	"x_0",		"mu02",	 	"thresh",	"C2", 	"mu202",	"g1",	"g2"};
double par_start[]	= {3.0e+01,	3.3e-01,	2.0e+00,	4.0,		5,		1.26,	2.0,		0.5,	0.5};
double par_error[]	= {	5.0,	0.05,		0.1,		1.0,		0.5,		0.01,	0.1,		0.1,	0.1 };

double   par_min[]	= {	0.0,	0.00,		0.0,		1.0,		0,		0.1,	1.0 ,		0.0,	0.0};
double   par_max[] 	= {	80.0,	1.00,		100.0,		10,		10,		10.0,	10.0,		2.0,	2.0};
//#endif
#elif (MODEL==1)||(MODEL==3)
std::string  par_name[]	= {
"sigma_0",	"A_g",	"lambda_g",	"mu02",		"C1", 	"mu102",	"thresh",		"C2",	"mu202",	"g1",	"g2"};
double par_start[]  	= {
33.0,		1.0,	0.2,		4.0,		0.2, 	1.5,   		1,		1.26,	2.0,		0.1, 	0.1};  
double par_error[]  	= {
 5.00,  	0.10,	0.2, 		1.0,		1.0,	0.1,		0.5,		0.1,	0.1,		0.01,	0.01};
double   par_min[] 	= {
 0.00,  	0.00,	-5.0,		1.0,		0.01, 	1.0,		0,		0.5,	1.0,		0.0,	0.0};
double   par_max[] 	= {
50.00,		20.00,	5.0, 		10.0,		10.0,	10.0,		10,		10,	10.0,		2.0,	2.0};
#endif

//int parameter(const std::vector<double>& par,double (& sigpar)[],double (& sudpar)[]){
int parameter(const std::vector<double>& par,double *sigpar,double *sudpar){
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
	sigpar[j++]=2.56819*(par[i++]);//mb to GeV
	sigpar[j++]=(par[i++]);
	sigpar[j++]=(par[i++])*1.0e-4;

#else
	sigpar[j++]=2.56819*(par[i++]);
	sigpar[j++]=(par[i++]);
	//sigpar[1]=pow(fabs(par[1])*1.0e-4,par[2]);
	sigpar[j++]=(par[i++]);
#endif

#if (MU02==0 && ALPHA_RUN==1)
	sigpar[j++]=par[i++];
#else 
//	++i;
#endif
//	printf("SIGMA: %.2e %.2e %.2e ",sigpar[0],sigpar[1],sigpar[2]);
#if(MODEL==1||MODEL==3)
	sigpar[j++]=(par[i++]);
	sigpar[j++]=(par[i++]);//sqrt(fabs(sigpar[3]/par[4]));//rmax^2= C/mu02
//	printf(" %.2e %.2e ",sigpar[3],sigpar[4]);
#endif
#if THRESHOLD==-1
	sigpar[j++]=(par[i++]);
#endif
////////////////////////////SUDPAR////////////////////////////////
#if (MODEL==0)
	sudpar[k++]=(par[i++]);
	sudpar[k++]=(par[i++]);
	//printf("\tSUDAKOV: %.2e %.2e ",sudpar[0],sudpar[1]);	
	#if (SUDAKOV==2)
		sudpar[k++]=(par[i++]);
		sudpar[k++]=(par[i++]);
		//printf("%.2e %.2e ",sudpar[2],sudpar[3]);
	#endif

///////////////////////////////////////////////////////
#elif (MODEL==1)
	#if INDEPENDENT_C==1
		sudpar[k++]=(par[i++]) ;
	#else 
		sudpar[k++]=(par[(i++)-2]);
	#endif

	#if INDEPENDENT_RMAX==1
		sudpar[k++]=(par[i++]);
	#else
		sudpar[k++]=(par[(i++)-2]);//mu02 is shared
	#endif
	//printf("%.2e %.2e ",sudpar[0],sudpar[1]);

//////////////////////////////////////////////////////
#if (SUDAKOV==2)
	sudpar[k++]=(par[i++]);
	sudpar[k++]=(par[i++]);
	//printf("%.2e %.2e ",sudpar[2],sudpar[3]);
#endif
#endif
//	printf("\n");
//	printf("sig %d  sud %d\n",j,k );
	//printf("\n");
	//if(i!=N_PAR){
	//	printf("parameter number mismatch counted %d, N_PAR= %d\n",i,N_PAR);
	//	exit(0);
	//}
	printf("sigpar: ");
	for(int l=0;l<j;l++){
		printf(" %.3e  ",sigpar[l]);
	}
	printf("\tsudpar: ");
	for(int l=0;l<k;l++){
		printf("%.3e  ",sudpar[l]);
	}printf("\n");
	return 0;
}
