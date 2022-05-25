#include<stdio.h>
#include<math.h>
#include<time.h>

#include"./control.h"
#include"./constants.h"


#include"simpson-integral.h"
#include"dipole-cross-section.h"
#include"photon-wave-function.h"

#include"DIS-cross-section.h"

#include"cfortran.h"
#include"../minuit.h"
#include"./read-and-fit.h"

//#define TEST 2

double arglist[10];
/* GBW Starting parameter values, errors: sigma_0,lambda,x_0, C, mu2,g1 */
double par_start[7]	= {	23.0,	0.29,	3.0e-4,		1.26,	4.0,	0.2,	0.8};
double par_error[7]	= {	1.0,	0.05,	0.1e-04,	0.01,	0.01,	0.01,	0.01 };

double   par_min[7]	= {	0.0,	0.00,	0.0,		0.01,	1.0 ,	0.0,	0.0};
double   par_max[7] 	= {	80.0,	1.00,	1.0,		20.0,	20.0,	2.0,	2.0};

int main(int argc, char* argv[]){
	int error_flag = 0;
	load_data();
	
#if TEST==1
	double resval=0.0;	
	double grad[7];
	fcn(7, grad, &resval , par_start,0, &dum_func );
	printf("fcn returns %f,\n",resval); 
	return(0);
#endif
	
	/* Initialize Minuit */
	MNINIT(5,6,7);
	/* Parameters definition */
#if ((MODEL==0)||( MODEL==2 ))  
	MNPARM(1,"sigma_0", par_start[0],par_error[0],par_min[0],par_max[0],error_flag);
	MNPARM(2," lambda", par_start[1],par_error[1],par_min[1],par_max[1],error_flag);
	MNPARM(3,"    x_0", par_start[2],par_error[2],par_min[2],par_max[2],error_flag);
	#if MODEL==2 
	MNPARM(4,"     C", par_start[3],par_error[3],par_min[3],par_max[3],error_flag);
    	MNPARM(5,"  bmax", par_start[4],par_error[4],par_min[4],par_max[4],error_flag);
	MNPARM(6,"    g1", par_start[5],par_error[5],par_min[5],par_max[5],error_flag);
	MNPARM(7,"    g2", par_start[6],par_error[6],par_min[6],par_max[6],error_flag);
	#endif
	
#elif MODEL==1
	MNPARM(1," sigma_0", par_start[0],par_error[0],par_min[0],par_max[0],error_flag);
	MNPARM(2,"     A_g", par_start[1],par_error[1],par_min[1],par_max[1],error_flag);
	MNPARM(3,"lambda_g", par_start[2],par_error[2],par_min[2],par_max[2],error_flag);
	MNPARM(4,"       C", par_start[3],par_error[3],par_min[3],par_max[3],error_flag);
	MNPARM(5,"     mu0", par_start[4],par_error[4],par_min[4],par_max[4],error_flag);
#endif

   

    arglist[0] = STRATEGY;

    /* Set strategy to STRategy from main.h, 0-fast, 1-default, 2-precise */
    MNEXCM(fcn,"SET STR",arglist, 1,error_flag,0);

    /* Minimalization procedure MIGRAD */
    MNEXCM(fcn,"MIGRAD",0,0,error_flag,0);
    
    return(0);
 
}
