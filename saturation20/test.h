/*  
 *
 */

//#include "functions.h"
#include "graph.h"
//#include "diffraction.h"
#include <stdlib.h>

/*******************************************************************************
* Check W2
*******************************************************************************/
void w2_test (void) {

    int i;
    double tx, tq2,tw,tdum;
    FILE  *datafile; 
    datafile = fopen("h1_epj_c21.dat","r");    
    for(i=0;i<13;i++) {
        fscanf(datafile,"%lE  %lE %lE %lE",&tq2,&tx,&tdum,&tdum ); 
        tw = sqrt(tq2*(1.0/tx-1.0));
        printf("%5.1f %f    %f\n",tq2,tx,tw); 
    }

    fclose(datafile);
}

/*******************************************************************************
* The function comparing my xgpdf with K.Golec-Biernat result
*******************************************************************************/
void comp_xgpdf_mg (void) {

    double xval[100];      /* Array for Bjoreken x values */
    double Qval;           /* Q^2 value */

    double q1valg[100];    /* K.G.-B. results at Q^=1    */
    double q10valg[100];   /* K.G.-B. results at Q^=10   */
    double q100valg[100];  /* K.G.-B. results at Q^=100  */
    double q1000valg[100]; /* K.G.-B. results at Q^=1000 */

    double kgb_value;
    double kgb_my_diff;
    double xgpdf_value;
    double ratio;

    //double xgpdfcheb; 
    int i;                /* Count */

    FILE  *datafile; 

    /* Open file with K.G-B. results */
    datafile = fopen("fig1a.dat","r");    

    /* Read the values of x and Q^2 together with K.G.-B. results */
    for(i=0;i<100;i++) {
        fscanf(datafile," %lE %lE %lE %lE %lE", 
               &xval[i], &q1valg[i], &q10valg[i], &q100valg[i], &q1000valg[i]);
    }

    fclose(datafile);

    /* Data from 'figa1.dat' were calculated for parameters below.
       If you want to change them do it at this place             */

    A_g  = 1.0;  
    beta = 6.0;
    Q0   = 1.0;
    n_f  = 3;
    Lambda2 = 0.0625;


    /* Print info */
    system("clear");   /* Clear sreen */
    printf("-----------------------------------------------------------------");
    printf("--------------\n");
    printf(" TESTING MODULE  -  xgpdf(x,Q^2) ");
    printf(" - exact value (not Chebyshev approximation)\n");
    printf("-----------------------------------------------------------------");
    printf("--------------\n");
    printf("\n");
    printf("Initial condition: xg(x,1) = %1.1f*x^%1.1f*(1-x)^%1.1f\n",
           A_g,fabs(-lambda),beta-1); 
    printf("\n");
    printf("Q_0^2  = %1.1f\n",Q0);
    printf("n_f = %d\n",n_f);
    printf("Lambda^2_QCD = %f\n",Lambda2);
    printf("\n");


    /* Ask for the value of Q^2 */
    printf("Choose the value of Q^2 from 1,10,100,1000\n");
    printf("Q^2 = ");
    scanf("%lf",&Qval);

 
    /* List comparison of my xpdf and K.-G.-B. xpdf */
    for(i=0;i<100;i++) {
        switch((int) Qval) {
            case 1:
            kgb_value = q1valg[i]; 
            xgpdf_value = xgpdf(xval[i],Qval);
            kgb_my_diff = fabs(xgpdf_value-kgb_value);
            break;
            case 10:
            kgb_value = q10valg[i]; 
            xgpdf_value = xgpdf(xval[i],Qval);
            kgb_my_diff = fabs(xgpdf_value-kgb_value);
            break;
            case 100:
            kgb_value = q100valg[i]; 
            xgpdf_value = xgpdf(xval[i],Qval);
            kgb_my_diff = fabs(xgpdf_value-kgb_value);
            break;
            case 1000:
            kgb_value = q1000valg[i]; 
            xgpdf_value = xgpdf(xval[i],Qval);
            kgb_my_diff = fabs(xgpdf_value-kgb_value);
            break;
            default:
            xgpdf_value = xgpdf(xval[i],Qval);
            break;

        }
       
        ratio = 200*kgb_my_diff/(kgb_value+xgpdf_value);
	printf(" %e %e %e  %.1e %c \n",xval[i],kgb_value,xgpdf_value,ratio,'%');
    }
}

/*******************************************************************************
* The function comparing exact xgpdf with Chebyshev approximation
*******************************************************************************/
void comp_xgpdf_ec (void) {

    FILE *out_file;
    //double Qval;           /* Q^2 value */


    double xgpdf_exact;
    double xgpdf_cheb;
    double exact_cheb_diff;

    int i;                /* Count */


    readdata();

    /* Data from 'figa1.dat' were calculated for parameters below.
       If you want to change them do it at this place             */

    A_g  = 1.0;  
    beta = 6.0;
    Q0   = 1.0;
    n_f  = 3;



    /* Print info */
    system("clear");   /* Clear sreen */
    printf("-----------------------------------------------------------------");
    printf("--------------\n");
    printf(" TESTING MODULE  -  xgpdf(x,Q^2) ");
    printf(" - exact vs. Chebyshev approximation\n");
    printf("-----------------------------------------------------------------");
    printf("--------------\n");
    printf("\n");
    printf("Initial condition: xg(x,1) = %1.1f*x^%1.1f*(1-x)^%1.1f\n",
           A_g,fabs(-lambda),beta-1); 
    printf("\n");
    printf("Q_0^2  = %1.1f\n",Q0);
    printf("n_f = %d\n",n_f);
    printf("\n");


    /* Perform chebyshev approximation for xgpdf */
    chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);
 

    /* List comparison of my xpdf and K.-G.-B. xpdf */
    out_file = fopen("exact_approx_cheb.txt","w");
    
    fprintf(out_file,"    xmod        xgpdf_exact  xgpdf_cheb   difference\n");
    fprintf(out_file,"\n");

    for(i=0;i<NDATA;i++) {


	xgpdf_exact = xgpdf(xmoddata[i],q2data[i]);
        xgpdf_cheb  = chebev(xmin,xmax,Qmin,Qmax,MX,MQ,coef,
                             xmoddata[i],q2data[i]);
        exact_cheb_diff = fabs(xgpdf_exact - xgpdf_cheb);
        
        fprintf(out_file," %e  %e %e  %e\n",xmoddata[i],
                xgpdf_exact,xgpdf_cheb,exact_cheb_diff);;

	printf(" %e  %e %e  %e\n",xmoddata[i],
               xgpdf_exact,xgpdf_cheb,exact_cheb_diff);;
    }

   fclose(out_file);
}

/*******************************************************************************
* The function comparing reduced cross sections with two methods of integration
*******************************************************************************/
void comp_int (void) {

    double sigma_dadmul;
    double sigma_simps2d;
    double chi_temp = 0.0;
    double chi2_dadmul = 0.0;
    double chi2_simps2d = 0.0;
    int i;                /* Counter */

    readdata();
    calc_ydata();

    /* Print info */
    system("clear");   /* Clear sreen */
    printf("-----------------------------------------------------------------");
    printf("--------------\n");
    printf(" TESTING MODULE  -  red_cs (x, Q^2) ");
    printf(" - dadmul vs. simps2d\n");
    printf("-----------------------------------------------------------------");
    printf("--------------\n");
    printf("\n");


    switch (model) {
	case 0:
	sigma_0  = par_start[0]; 
	lambda   = par_start[1]; 
	x_0      = par_start[2]; 
	break;

	case 1:
	sigma_0  = bgk_parst[0]; 
	A_g      = bgk_parst[1]; 
	lambda_g = bgk_parst[2]; 
	C        = bgk_parst[3]; 
	mu02     = bgk_parst[4]; 
	break;
    }


    /* Perform chebyshev approximation for xgpdf */
    if (model == 1)
    chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);

    FILE *out_file;
    out_file = fopen("result__.txt","w");
    /* List comparison of */
    for(i=0;i<NDATA;i++) {

         uif_int = 0;         
         sigma_dadmul = sigma(xmoddata[i],q2data[i],ydata[i],bgk_parst);
         chi_temp = (sigma_dadmul - csdata[i])/erdata[i];
         chi2_dadmul  += chi_temp*chi_temp;

         uif_int = 1;         
         //sigma_simps2d = sigma(xmoddata[i],q2data[i],ydata[i],par_start);
         sigma_simps2d = sigma(xmoddata[i],q2data[i],ydata[i],bgk_parst);
         chi_temp = (sigma_simps2d- csdata[i])/erdata[i];
         chi2_simps2d += chi_temp*chi_temp;

         printf("%.1e  %e  %e  %e  %e\n",q2data[i],xmoddata[i],csdata[i],
                sigma_dadmul, sigma_simps2d);
         fprintf(out_file,"%.1e  %e  %e  %e  %e\n",
                q2data[i],xmoddata[i],csdata[i],
                sigma_dadmul, sigma_simps2d);

    }

    fclose(out_file);
    printf("---------------------------------------------------------\n");
    printf("chi^2_dadmul  =  %e\n",chi2_dadmul);
    printf("chi^2_simps2d =  %e\n",chi2_simps2d);
}


/*******************************************************************************
* The function  testing dadmul integration
*******************************************************************************/
void dadmul_test (void) {

    xmod = 1.0e-03;
    Q2 = 10;
    y = 0.07;

    /* Integration limits */
    double Rmax_test = 1.0e+03;
    double A[2] = {1.0e-10,1.0e-10}; /* Lower limits of integration  {r,z} */ 
    double D[3] = {10.0,1.0,100};
    double B[2] = {Rmax_test,0.5};   /* Upper limits of integration  {r,z} */ 


    /* dadmul integration variables */
    int    dim = 2;            /* Dimension of integral */
    int    minpts = 1.7e+01;
    int    maxpts = 1.0e+05;
    //int    maxpts = 1.0e+04;
    double eps = 1.0e-06;     /* Relative accuracy */
    //int    iwk = 2100;
    int    iwk = 21000;
    double wk[iwk];
    double result;             /* The result of integration */
    double relerr;             /* The relative accuracy of result */ 
    int    nfnevl;             /* The number of function eval. performed */
    int    ifail;              /* If '0' normal exit relerr < eps */
    double temp = 0.0,temp1;

    /* Set model parameters according to sigma() argument 'par' */    
    sigma_0  = bgk_parst[0]; 
    A_g      = bgk_parst[1];
    lambda_g = bgk_parst[2];
    C        = bgk_parst[3]; 
    mu02     = bgk_parst[4];

    chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);

    //B[0] = 1.0;
    B[0] = D[0];
    dadmul_(&uif_bgk, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
	    &result, &relerr, &nfnevl, &ifail);
    temp += result;
/*
    //A[0] = 1.0;
    //B[0] = 10.0;
    A[0] = D[0];
    B[0] = D[1];

    dadmul_(&uif_bgk, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
	    &result, &relerr, &nfnevl, &ifail);
    temp += result;

    A[0] = D[1];
    B[0] = D[2];

    dadmul_(&uif_bgk, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
	    &result, &relerr, &nfnevl, &ifail);
    temp += result;
*/
    A[0] = D[0];
    B[0] = Rmax_test;

    dadmul_(&uif_bgk, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
	    &result, &relerr, &nfnevl, &ifail);
    temp += result;


    A[0] = 1.0e-10;
    A[1] = 1.0e-10;
    B[0] = Rmax_test;

    temp1 = simps2d(A,B,1.0e-05, uif_bgk);

    printf("dadmul     simps2d\n"); 
    printf("%e       %e\n", temp, temp1); 

}

/*******************************************************************************
* The function  testing dadmul integration
*******************************************************************************/
void psi_charm_test(void) {


    FILE *out_file;
    double chi2_lc;
    double r_i;
    double z_i;
    int i,j;


    /* Print info */
    system("clear");   /* Clear sreen */
    printf("-----------------------------------------------------------------");
    printf("--------------\n");
    printf(" TESTING MODULE  -  red_cs (x, Q^2) ");
    printf(" - dadmul vs. simps2d\n");
    printf("-----------------------------------------------------------------");
    printf("--------------\n");
    printf("\n");

    out_file = fopen("result_psi.txt","w");

    /* List comparison of */

    for(i=0;i<5;i++) {

        r_i = pow(10,i);

	for(j=1;j<5;j++) {
        z_i = 0.1*j;

         chi2_lc = psisqlc(r_i,z_i);
         printf("%e  %e  %e\n",r_i,z_i,chi2_lc);

/*
         fprintf(out_file,"%.1e  %e  %e  %e  %e\n",
                q2data[i],xmoddata[i],csdata[i],
                sigma_dadmul, sigma_simps2d);
*/
         }
    }

    fclose(out_file);

}

