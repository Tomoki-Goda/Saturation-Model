//#include "diffraction.h"
#include "functions.h"
#include <math.h>

double arg1,arg2;
//double arg1,arg2;
//double drmin = 8.0;
double drmin = 1.0e-04;
//double drmax = 8.0;
double drmax = 1.0e+04;
double dumin = 1.0e-04;
double dumax = 1.0e+04;

double dzmin;
double dzmax = 0.5;
double deps  = 1.0e-05;

double bessk(int n, double x);
double bessj(int n, double x);
double m2;


/*******************************************************************************
* The function plotting 
*******************************************************************************/
void graph_diff_xph1zeus(void) {

    //FILE *fout; 
    // int i,k,j;

    //double q2diff[12]={3.5,5.0,6.5,8.5,12.0,15.0,20.0,25.0,35.0,45.0,60.0,90.0};
    //double dbeta_graph[7]= {0.01,0.04,0.1,0.2,0.4,0.65,0.9};

    //double dxp,FT,FL,Fg = 0.0;

    //                     /*{0.949,0.750,0.4290,0.1820,0.0632,0.0291},\ */
    //double npoints; /* Number of points used to plot */

    //if (flavour == 1)
    //    chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);


    //npoints = 8;

    //fout = fopen("dxp_model_h1zeus.dat","w");    
    ////fout = fopen("diff_xpMODEL.dat","w");    

    //for(i=0;i<12;i++) {        /* Row , 12 */
    //    for(k=0;k<7;k++) {     /* Column, 7 */
    //        for(j=0;j<npoints;j++) {
    //    	//dxp = 0.01;
    //    	dxp = 0.05*pow(10,-0.3*j);
    //    	xmod = dxp;
    //    	FT = FD_T(q2diff[i],dbeta_graph[k],xp,m_fsq);
    //    	FL = FD_L(q2diff[i],dbeta_graph[k],xp,m_fsq);
    //    	Fg = FD_gluon(q2diff[i],dbeta_graph[k],xp);

    //    	//printf("%e  %f\n",dxp, dbeta_graph[i][k]);
    //    	printf("%4.1f  %4.2f     %e  %e  %e\n",
    //                    q2diff[i], dbeta_graph[k],dxp, FT+FL,Fg);
    //    	fprintf(fout,"%e  %e  %e\n",dxp, FT+FL,Fg);
    //        }
    //        fprintf(fout,"\n"); 
    //    }
    //    //printf("Row %d finished \n",i+1); 
    //    printf("Column %d finished\n",k+1); 
    //}
    //fclose(fout);
}


/*******************************************************************************
* Plot of alpha_s * gluon
*******************************************************************************/
void graph_alphagluon () {
   
    int i;
    double mu2,r;
    double s_gbw1,s_bgk1,s_gbw2,s_bgk2;
    FILE *out_file;

    switch (model) {
    case 0:
    break;
    case 1:
	sigma_0  = bgk_parst[0]; 
	A_g      = bgk_parst[1]; 
	lambda_g = bgk_parst[2]; 
	C        = bgk_parst[3]; 
	mu02     = bgk_parst[4]; 
    break;
    }

    


    out_file = fopen("alpha_gluon.dat","w");
 
    for (i=0;i<54;i++) {
        r = 0.00001*pow(10,0.1*i);
	mu2 = C/(r*r) + mu02;

	xmod = 0.0001;
	s_bgk1 = 0.389379*pi*pi/3.0*alpha_s(mu2)*xgpdf(xmod,mu2);
	s_gbw1 = 18.8*0.25*pow(xmod/0.0003,-0.29);

	xmod = 0.01;
	s_bgk2 = 0.389379*pi*pi/3.0*alpha_s(mu2)*xgpdf(xmod,mu2);
	s_gbw2 = 18.8*0.25*pow(xmod/0.0003,-0.29);

        fprintf (out_file,"%e   %e   %e   %e   %e\n",
                 r, s_bgk1,s_bgk2, s_gbw1,s_gbw2);
    }

    fclose(out_file);
   
}

/*******************************************************************************
* Plot of alpha_s * gluon
*******************************************************************************/
void graph_gluon () {
   
    int i,k, npoints;
    double s_bgk[4], x_Bj;
    double mu2[4] = {2.0,5.0,20.0,100.0};
    FILE *out_file;

    model = 1;

    sigma_0  = bgk_parst[0]; 
    A_g      = bgk_parst[1]; 
    lambda_g = bgk_parst[2]; 
    C        = bgk_parst[3]; 
    mu02     = bgk_parst[4]; 

    out_file = fopen("gluon.dat","w");
 
    npoints = 30;
    for (i=0;i<npoints;i++) {
	x_Bj = pow(10,-(4.0*i/(npoints-1)+1.0)); 

        for(k=0;k<4;k++) 
	    s_bgk[k] = xgpdf(x_Bj,mu2[k]);

        printf ("%.2e   %e   %e   %e   %e\n",
                 x_Bj, s_bgk[0],s_bgk[1],s_bgk[2],s_bgk[3]); 
        fprintf (out_file,"%e   %e   %e   %e   %e\n",
                 x_Bj, s_bgk[0],s_bgk[1],s_bgk[2],s_bgk[3]); 
    }
    fclose(out_file);
}


/*******************************************************************************
* Plot of F_L
*******************************************************************************/
void graph_fl () {
   
    int i,k,npoints;
    double Q2_g,x_g;
    FILE *out_file;

    double q2_fl[5] = {0.75,1.35,2.2,4.2,7.5};
    switch (model) {
    case 0:
    break;
    case 1:
	sigma_0  = bgk_parst[0]; 
	A_g      = bgk_parst[1]; 
	lambda_g = bgk_parst[2]; 
	C        = bgk_parst[3]; 
	mu02     = bgk_parst[4]; 
	chebft(xmin,xmax,Qmin,Qmax,MX,MQ,*coef, &xgpdf);
    break;
    }

    


    out_file = fopen("fl_bgk.dat","w");
 
    npoints = 10;

    dataform = 2;
    printf("dataform = %d\n",dataform);
    printf("flavour  = %d\n",flavour);

    for (k=0;k<5;k++) {
	Q2_g = q2_fl[k];
	for (i=0;i<npoints;i++) {
	    x_g = 0.1*pow(10,-(i*0.4+1.38));
	    printf("%e  %e\n", x_g,sigma(x_g,Q2_g,0,bgk_parst));
	    fprintf(out_file, "%e  %e  %e\n", x_g,
				  sigma(x_g,Q2_g,0,bgk_parst),
				  sigma_c(x_g,Q2_g,0,bgk_parst));
	}
	fprintf(out_file,"\n");
    }

    fclose(out_file);
   
}

/*******************************************************************************
* Plot of F_L as a function of Q2
*******************************************************************************/
void graph_flq2 (int model_graph,int flavour_graph) {
   
    int k,npoints;
    double Q2_g,x_g;
    double W = 276;
    FILE *out_file;
    double par_temp[5];

    model = model_graph;
    flavour = flavour_graph;

    switch (model) {
    case 0:
	model = 0;
        xbj_mod = 0; 
        switch (flavour) {
	    case 1:
	    par_temp[0] = 29.12; 
	    par_temp[1] = 0.277; 
	    par_temp[2] = 0.41e-04; 

	    printf("Calculating F_L(Q^2) in GBW model");
	    printf(" and saving it in flq2_gbwc.dat\n");

	    out_file = fopen("flq2_gbwc.dat","w");    
	    break;
	    case 0:
	    par_temp[0] = 23.03; 
	    par_temp[1] = 0.288; 
	    par_temp[2] = 3.04e-04; 

	    printf("Calculating F_L(Q^2) in GBW model");
	    printf(" and saving it in flq2_gbwl.dat\n");

	    out_file = fopen("flq2_gbwl.dat","w");    
	    break;
        }
    break;
    case 1:
        model = 1;
	sigma_0  = bgk_parst[0]; 
	A_g      = bgk_parst[1]; 
	lambda_g = bgk_parst[2]; 
	C        = bgk_parst[3]; 
	mu02     = bgk_parst[4]; 

	par_temp[0] = sigma_0; 
	par_temp[1] = A_g; 
	par_temp[2] = lambda_g; 
	par_temp[3] = C ; 
	par_temp[4] = mu02 ; 

        printf("\nCalculating F_L(Q^2) in BGK model");
        printf(" and saving it in flq2_bgk.dat\n");

	out_file = fopen("flq2_bgk.dat","w");
	chebft(xmin,xmax,Qmin,Qmax,MX,MQ,*coef, &xgpdf);
    break;
    }
 
    npoints = 33;

    dataform = 2;
    printf("dataform = %d (F_L = 2)\n\n",dataform);

    for (k=0;k<npoints;k++) {
        Q2_g = pow(10,0.1*k-0.3);
	x_g = Q2_g/(W*W+Q2_g);
	    printf("%e  %e\n", Q2_g,sigma(x_g,Q2_g,0,par_temp));
	    fprintf(out_file, "%e  %e  %e\n", 
                    Q2_g,sigma(x_g,Q2_g,0,par_temp),
                         sigma_c(x_g,Q2_g,0,par_temp)
                         + sigma_b(x_g,Q2_g,0,par_temp));
    }
    fclose(out_file);
}

/*******************************************************************************
* Critical line
*******************************************************************************/
void graph_crit_line(int model_graph) {
   
    int k,npoints;
    double Q2_bgk,Q2_gbw,logx;
    FILE *out_file;

    model = model_graph;

    sigma_0  = bgk_parst[0]; 
    A_g      = bgk_parst[1]; 
    lambda_g = bgk_parst[2]; 
    C        = bgk_parst[3]; 
    mu02     = bgk_parst[4]; 

    printf("BGK\n");
    //out_file = fopen("crit_linel.dat","w");
    out_file = fopen("crit_linecb.dat","w");
    chebft(xmin,xmax,Qmin,Qmax,MX,MQ,*coef, &xgpdf);
 
    npoints = 25;


    for (k=0;k<npoints;k++) {
        logx = -0.18*k-2.0;
        //logx = -0.18*k-2.5;
	xmod = pow(10,logx);
        Q2_bgk = find_zero (0,50, &crit_func);  /* BGK light */
        //Q2_gbw = pow(0.41e-04/xmod,0.277);    /* GBW charm */
        Q2_gbw = pow(0.1632e-04/xmod,0.2197);      /* Gregory */
        //Q2_gbw = pow(3.04e-04/xmod,0.288);      /* GBW light */
        //Q2_gbw= find_zero (0,50, &crit_func_frozen);  /* BGK light */
	    printf("%e  %e  %e\n",Q2_bgk,Q2_gbw,-logx);
	    fprintf(out_file,"%e  %e  %e\n",Q2_bgk,Q2_gbw,-logx);
    }
    fclose(out_file);
}

/*******************************************************************************
* Saturation scale
*******************************************************************************/
void graph_sat_scale(void) {
   
    int k,npoints;
    double Q2_gbw_3,Q2_gbw_4,logx;
    double Q2_bgk_full,Q2_bgk_frozen;
    FILE *out_file;


    sigma_0  = bgk_parst[0]; 
    A_g      = bgk_parst[1]; 
    lambda_g = bgk_parst[2]; 
    C        = bgk_parst[3]; 
    mu02     = bgk_parst[4]; 

    printf("BGK\n");
    out_file = fopen("sat_scale.dat","w");
    chebft(xmin,xmax,Qmin,Qmax,MX,MQ,*coef, &xgpdf);
 
    npoints = 20;


    for (k=0;k<npoints;k++) {
        logx = -0.18*k-2.0;
	xmod = pow(10,logx);
        Q2_gbw_3      = pow(3.04e-04/xmod,0.288);             /* GBW light  */
        Q2_gbw_4      = pow(0.41e-04/xmod,0.277);             /* GBW charm  */
        Q2_bgk_full   = find_zero (0,50, &crit_func);         /* BGK full   */
        Q2_bgk_frozen = find_zero (0,50, &crit_func_frozen);  /* BGK frozen */

	   printf("%e  %e  %e %e  %e\n",
                   xmod,Q2_gbw_3,Q2_gbw_4,Q2_bgk_full,Q2_bgk_frozen);
	   fprintf(out_file,"%e  %e  %e %e  %e\n",
                   xmod,Q2_gbw_3,Q2_gbw_4,Q2_bgk_full,Q2_bgk_frozen);
    }
    fclose(out_file);
}

/*******************************************************************************
* Saturation scale nf= 0, nf = 5
*******************************************************************************/
void graph_sat_scale_05(void) {
   
    int k,npoints;
    double logx;
    double Q2_bgk_full,Q2_bgk_frozen;
    FILE *out_file;


    sigma_0  = bgk_parst[0]; 
    A_g      = bgk_parst[1]; 
    lambda_g = bgk_parst[2]; 
    C        = bgk_parst[3]; 
    mu02     = bgk_parst[4]; 

    printf("BGK\n");
    out_file = fopen("sat_scale0.dat","w");
    chebft(xmin,xmax,Qmin,Qmax,MX,MQ,*coef, &xgpdf);
 
    npoints = 20;


    for (k=0;k<npoints;k++) {
        logx = -0.18*k-2.0;
	xmod = pow(10,logx);
        Q2_bgk_full   = find_zero (0,50, &crit_func);         /* BGK full   */
        Q2_bgk_frozen = find_zero (0,50, &crit_func_frozen);  /* BGK frozen */

	   printf(" %e %e  %e\n",
                   xmod,Q2_bgk_full,Q2_bgk_frozen);
	   fprintf(out_file,"%e %e  %e\n",
                   xmod,Q2_bgk_full,Q2_bgk_frozen);
    }
    fclose(out_file);
}


/*******************************************************************************
* Unintegrated gluon distribution graph
*******************************************************************************/
void alpha_fpdf_graph(void) {
   
    #define X_N 3        /* Number of lines */
    #define NPOINTS 60   /* Number of points used to plot */

    int i,k;
    double x_g;
    double a_fpdf_0[X_N][NPOINTS],a_fpdf_1[X_N][NPOINTS];
    double l2_g[NPOINTS];
    FILE *out_file;

    /******* GBW model calculations ********/

    model = 0;  /* Information for alpha_fpdf(x,l2) */
	sigma_0 = 29.12;  /* GBW model with charm parameters */
        lambda  = 0.277;
        x_0     = 0.41e-04;  

	printf("\nCalculating alpha_fpdf(x,l^2) in GBW model...");

    for (i=0;i<NPOINTS;i++) {
        l2_g[i] = 1.0e-02*pow(10,0.1*i);

	for (k=0;k<X_N;k++) {
	    xmod = 1.0e-04*pow(10,k);
            x_g = xmod;
	    a_fpdf_0[k][i] = alpha_fpdf(x_g,l2_g[i]);
        }
    }

    /******* BGK model calculations ********/

    model = 1; /* Information for alpha_fpdf(x,l2) */
    sigma_0  = bgk_parst[0]; 
    A_g      = bgk_parst[1]; 
    lambda_g = bgk_parst[2]; 
    C        = bgk_parst[3]; 
    mu02     = bgk_parst[4]; 

    printf("\nCalculating alpha_fpdf(x,l^2) in BGK model...\n");

    chebft(xmin,xmax,Qmin,Qmax,MX,MQ,*coef, &xgpdf);
    for (i=0;i<NPOINTS;i++) {
        l2_g[i] = 1.0e-02*pow(10,0.1*i);

	for (k=0;k<X_N;k++) {
	    xmod = 1.0e-04*pow(10,k);
            x_g = xmod;
	    a_fpdf_1[k][i] = alpha_fpdf(x_g,l2_g[i]);
        }
    }

    /******* BGK massive original model calculations ********/
/*
    model = 1; 
    sigma_0  = 23.0; 
    A_g      = 1.2; 
    lambda_g = 0.28; 
    C        = 0.26; 
    mu02     = 0.52; 

    n_f = 3;
    printf("\nCalculating alpha_fpdf(x,l^2) in BGK model...\n");

    chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);
    for (i=0;i<NPOINTS;i++) {
        l2_g[i] = 1.0e-02*pow(10,0.1*i);

	for (k=0;k<X_N;k++) {
	    xmod = 1.0e-04*pow(10,k);
            x_g = xmod;
	    a_fpdf_0[k][i] = alpha_fpdf(x_g,l2_g[i]);
        }
    }
*/
    /* Print the results on the screen and to the file */

    printf("Results saved in 'a_fpdf.dat'. \n \n");
        //printf("l^2           x=10^{-4}      x=10^{-3}      x=10^{-2}\n\n");

    out_file = fopen("a_fpdf.dat","w");

    for (i=0;i<NPOINTS;i++) {
	fprintf(out_file,"%.3e  % e  % e  % e  % e  % e  % e\n",
             l2_g[i], a_fpdf_0[0][i],a_fpdf_0[1][i],a_fpdf_0[2][i],
                   a_fpdf_1[0][i],a_fpdf_1[1][i],a_fpdf_1[2][i]);
	//printf("%.3e  % e  % e  % e\n",
        //     l2_g[i], a_fpdf_0[0][i],a_fpdf_0[1][i],a_fpdf_0[2][i]);
    }
    fclose(out_file);
}

/*******************************************************************************
* The function plotting F_2_charm
*******************************************************************************/
void graph_dipole_scaling (void) {

    FILE *fout; 
     int i,k;

    double npoints; /* Number of points used to plot */
    double r,tau,rmin = 0.001;
    double delta;
    double Q2_scal;
    double dcsr[5];
    //double dcsr[5],dcsx[5];

    sigma_0  = bgk_parst[0]; 
    A_g      = bgk_parst[1]; 
    lambda_g = bgk_parst[2]; 
    C        = bgk_parst[3]; 
    mu02     = bgk_parst[4]; 

    npoints = 70;

    chebft(xmin,xmax,Qmin,Qmax,MX,MQ,*coef, &xgpdf);
    delta = (log10(30)-log10(rmin))/npoints; /* Logarithmic x-scale */

	fout = fopen("dipol_scaling.dat","w");    
	for(i=0;i<=npoints;i++) {
            tau = rmin*pow(10,delta*i);
            for(k=0;k<5;k++) {
                xmod = pow(10,-(k+2)); /* x = 10^-2...10^-6 */
                    Q2_scal = find_zero (0,50, &crit_func_frozen); 
                    r = sqrt(tau/Q2_scal); 
		    dcsr[k] = sigma_bgk(r);
		    //dcsr[k] = sigma_scaling(r,Q2_scal);
            }

	    fprintf(fout,"%e %e %e %e %e %e\n",
                   tau, dcsr[0], dcsr[1], dcsr[2], dcsr[3], dcsr[4]);
        }
	fclose(fout);

}
