/*  
 *
 */

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
* 
*******************************************************************************/
double phi_0_integrand (double *r)  {   

    double z = arg1;
    double dipol_cs;
    double value;

    double  Qsq_bar =  z*(1-z)*Q2+m2;
    double  Qsq2    =  sqrt(Qsq_bar)*r[0];

    double  M2 = Q2*(1.0/dbeta-1.0);
    double  k2 = fabs(z*(1-z)*M2-m2);
    double  k  = sqrt(k2)*r[0];

    double bessel_k0 = dbesk0_(&Qsq2);
    double bessel_j0 = dbesj0_(&k);
    //printf("%f\n",Qsq_bar);
    //printf("%f\n",k);

    //xmod = xp;
    switch (model) {
        case 0:
	dipol_cs = sigma_gbw(r[0]); 
        break;
        case 1:
	//dipol_cs = sigma_bgk(r[0]); 
	dipol_cs = sigma_bgk_cheb(r[0]); 
        break;
    }
    if (dipol_cs<0)
	printf("Warning ! Negative gluon !\n");

    //printf("%f   %f\n",r[0],dipol_cs);

    value = r[0]*bessel_k0*bessel_j0*dipol_cs; 
    return value;
}

/*******************************************************************************
* 
*******************************************************************************/
double phi_1_integrand (double *r)  {   

    double z = arg2;
    double dipol_cs;
    double value;

    double Qsq_bar =  z*(1-z)*Q2+m2;
    double Qsq2    =  sqrt(Qsq_bar)*r[0];

    double M2 = Q2*(1.0/dbeta-1.0);
    double k2 = fabs(z*(1-z)*M2-m2);
    double k  = sqrt(k2)*r[0];

    //double bessel_k1 = bessk1(Qsq2);
    //double bessel_j1 = bessj1(k);
    double bessel_k1 = dbesk1_(&Qsq2);
    double bessel_j1 = dbesj1_(&k);

    //printf("%f\nz);
    //printf("%f  %f\n",z, dzmin);

    //xmod = xp;
    switch (model) {
        case 0:
	dipol_cs = sigma_gbw(r[0]); 
        break;
        case 1:
	//dipol_cs = sigma_bgk(r[0]); 
	dipol_cs = sigma_bgk_cheb(r[0]); 
        break;
    }
    if (dipol_cs<0)
	printf("Warning ! Negative gluon !\n");

    //printf("%f\n",dipol_cs);
    value = r[0]*bessel_k1*bessel_j1*dipol_cs; 
    return value;
}

/*******************************************************************************
* 
*******************************************************************************/
double phi_0 (double *z)  {   

    //extern double dgquad_();
    extern double simps_();
    //int N = 96;
    double dum1,dum2,dum3;
    double result;
   
    arg1 = z[0];
    simps_(&drmin,&drmax,&deps,&deps,&deps,
           &phi_0_integrand,&dum1,&result,&dum2,&dum3);

    //value = dgquad_(&phi_0_integrand,&drmin,&drmax,&N);
	
    return result;
}

/*******************************************************************************
* 
*******************************************************************************/
double phi_1 (double *z)  {   

    //extern double dgquad_();
    extern double simps_();
    //int N = 96;
    double dum1,dum2,dum3;
    double result;
   
    arg2 = z[0];
    simps_(&drmin,&drmax,&deps,&deps,&deps,
           &phi_1_integrand,&dum1,&result,&dum2,&dum3);

    //value = dgquad_(&phi_0_integrand,&drmin,&drmax,&N);
	
    return result;
}


/*******************************************************************************
*									       **                                    G L U O N S                               *
*									       *
*******************************************************************************/

/*******************************************************************************
* 
*******************************************************************************/
double phi_2_integrand (double *u)  {   

    double k = arg1;
    double z = arg2;
    double dipol_cs;
    double bessel_k2 = bessk(2,sqrt(z/(1.0-z))*u[0]);
    double bessel_j2 = bessj(2,u[0]);

    double value;

   //printf("phi_2_int: %f\n",arg1);
    switch (model) {
        case 0:
	dipol_cs = sigma_gbw(u[0]/k); 
        break;
        case 1:
	dipol_cs = sigma_bgk_cheb(u[0]/k); 
	//dipol_cs = sigma_bgk_cheb_mod(u[0]/k); 
        break;
    }

    if (dipol_cs<0)
	printf("Warning ! Negative gluon !\n");

    //printf("%f\n",dipol_cs);
    value = u[0]*bessel_k2*bessel_j2*dipol_cs; 
    return value;
}


/*******************************************************************************
* 
*******************************************************************************/
double phi_2 (double *k)  {   

    extern double simps_();
    double dum1,dum2,dum3;
    double z;
    double int_result;
    double value;
    
    arg1 = k[0];
    z = arg2;

    //printf("phi_2: %f\n",arg1);
    simps_(&dumin,&dumax,&deps,&deps,&deps,
          &phi_2_integrand,&dum1,&int_result,&dum2,&dum3);
 
    value = 2*k[0]*log((1.0-z)*Q2/(k[0]*k[0]))*int_result*int_result;
    return value;
}

/*******************************************************************************
* 
*******************************************************************************/
double FD_gluon_integrand (double *z)  {   

    extern double dadapt_();

    double factor1,factor2;
    int nseg = 1;
    double kmin,kmax;
    double reltol,abstol,result,reserr;
    double value;
    arg2 = z[0];

    kmin = 0.0;
    kmax = sqrt((1-z[0])*Q2);
    reltol = 1.0e-04;
    abstol = 0.0;

    factor1 = 1.0/((1.0-z[0])*(1.0-z[0])*(1.0-z[0]));
    factor2 = (1.0-dbeta/z[0])*(1.0-dbeta/z[0])+dbeta/z[0]*dbeta/z[0];
    dadapt_(&phi_2,&kmin,&kmax,&nseg,&reltol,&abstol,&result,&reserr);


    value = factor1*factor2*result;
    return value;
}

/*******************************************************************************
* 
*******************************************************************************/
double FD_gluon (double Q2_g, double dbeta_g, double xp_g)  {   

    extern double dgquad_();
    //int N = 64;
    int N = 96;

    double dalpha_s = 0.2;
    double ef_sum, norm_factor;
    double value;
    
    Q2 = Q2_g;
    dbeta = dbeta_g;
    //xp = Q2_g*(1.0/xp_g-1.0);
    xp = xp_g;
    //xmod = xp;

    norm_factor = 81*dbeta_g*dalpha_s/(hc2*hc2*512*pow(PI,5)*B_D);
    ef_sum = 2.0/3.0;
    dzmin = dbeta_g;
    dzmax = 1;

    //printf("%f\n",dzmin);
    //printf("%f\n",dzmax);
    value = norm_factor*ef_sum*dgquad_(&FD_gluon_integrand,&dzmin,&dzmax,&N);
    return value;
}
/*******************************************************************************
* 
*******************************************************************************/
double FD_L_integrand (double *z)  {   

    double value;
    
    value = z[0]*z[0]*z[0]*(1-z[0])*(1-z[0])*(1-z[0])*phi_0(z)*phi_0(z);
    return value;
}
/*******************************************************************************
* 
*******************************************************************************/
double FD_L (double Q2_L, double dbeta_L, double xp_L, double m2_L)  {   

    extern double dgquad_();
    //int N = 64;
    int N = 96;

    double M2_L;
    double ef_sum, norm_factor;
    double value;
    
    Q2 = Q2_L;
    dbeta = dbeta_L;
    xp = xp_L;
    m2 = m2_L;

    norm_factor = 3*pow(Q2_L,3)/(hc2*hc2*32*pow(PI,4)*dbeta_L*B_D);
    ef_sum = 2.0/3.0;
    M2_L = Q2*(1.0/dbeta-1.0);
    dzmin = 0.5*(1-sqrt(1-4.0*m2/M2_L));
    dzmax = 0.5;


    value = norm_factor*ef_sum*2*dgquad_(&FD_L_integrand,&dzmin,&dzmax,&N);
    return value;
}

/*******************************************************************************
* 
*******************************************************************************/
double FD_T_integrand (double *z)  {   

    double epsilon2 =  z[0]*(1-z[0])*Q2+m2;
    double value;
    value = z[0]*(1-z[0])
            *(epsilon2*(z[0]*z[0]+(1-z[0])*(1-z[0]))*phi_1(z)*phi_1(z)
                + m2*phi_0(z)*phi_0(z));
    //printf("%f\n",phi_0(z)*phi_0(z));
    //printf("%f\n",phi_1(z)*phi_1(z));
    return value;
}

/*******************************************************************************
* 
*******************************************************************************/
double FD_T (double Q2_T, double dbeta_T, double xp_T, double m2_T)  {   

    extern double dgquad_();
    //int N = 64;
    int N = 96;

    double  M2_T;
    //double temp;
    double ef_sum, norm_factor;
    double value;
    
    Q2 = Q2_T;
    dbeta = dbeta_T;
    xp = xp_T;
    //xmod = xp;
    m2 = m2_T;

    norm_factor = 3*pow(Q2_T,2)/(hc2*hc2*128*pow(PI,4)*dbeta_T*B_D);
    ef_sum = 2.0/3.0;

    M2_T= Q2*(1.0/dbeta-1);
    //temp = Q2*(1.0/dbeta-1);
    dzmin = 0.5*(1-sqrt(1-4.0*m2/M2_T));
    dzmax = 0.5;

    //double temp = Q2*(1.0/dbeta-1);
    //printf("%e   %e\n",dzmin,pow(Q2_T,2));
    //printf("%e   %e   %e   %e\n",m2,M2_T,dzmin,dzmin*(1-dzmin)*M2_T-m2);
/* NOWE POCZATEK 
    arg1 = 0.2;
    double xyz = 0.001;
    value = phi_0_integrand(&xyz);
    //printf("%f\n", xp);
 NOWE KONIEC */

/* TO JEST POPRAWNE */
    value = norm_factor*ef_sum*2*dgquad_(&FD_T_integrand,&dzmin,&dzmax,&N);
    return value;
}

/*******************************************************************************
* 
*******************************************************************************/
double sigma_DT_integrand (double *M2_T){

    double xpl,dbl,Q2l,m2l;
    
    xpl = (Q2+M2_T[0])/W2;
    dbl = Q2/(Q2+M2_T[0]);
    Q2l = Q2;
    m2l = m2;

    return FD_T (Q2l,  dbl, xpl, m2l);

}

double sigma_DT (double W2_DT, double Q2_DT, double Mmin, double Mmax, 
                 double m2_DT){

    double value;
    double deps;
    extern double dgauss_();
    Q2 = Q2_DT;
    W2 = W2_DT;
    m2 = m2_DT;
   
    deps = 100000;
    value = dgauss_(&sigma_DT_integrand,&Mmin,&Mmax,&deps);
    return value;
}

/*******************************************************************************
*									       **                                    G R A P H S                               *
*									       *
*******************************************************************************/

/*******************************************************************************
* The function plotting 
*******************************************************************************/
void graph_diff (void) {

    FILE *fout; 
     int i,k;

    double q2diff[4] = {2.7,6.0,14,27.0};
    double dbeta_graph = 0.0, FT,FL,Fg = 0.0;

    double npoints; /* Number of points used to plot */

    if (model == 1)
	chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);

    xp = 0.01;
    xmod = xp;

    npoints = 10;

    //printf("%e\n",FD_T(q2diff[0],0.6,xp,m_fsq));

    fout = fopen("graph_diff.dat","w");    

    for(i=0;i<=3;i++) {
	//Q2 = q2diff[i];
	for(k=0;k<npoints;k++) {
	    dbeta_graph = 0.105*k+0.01; /* x = 10^-1...10^-6 */
	    //x_Bj = pow(10,-(0.5*k+1)); /* x = 10^-1...10^-6 */

	    FT = FD_T(q2diff[i],dbeta_graph,xp,m_fsq);
                //+2.0/3.0*FD_T(Q2,dbeta,xp,m_ch);
	    FL = FD_L(q2diff[i],dbeta_graph,xp,m_fsq);
                //+2.0/3.0*FD_L(Q2,dbeta,xp,m_ch);
	    //printf("%e %e %e %e\n", 
	    Fg = FD_gluon(q2diff[i],dbeta_graph,xp);

	    fprintf(fout,"%e  %e  %e  %e  %e\n", 
			  dbeta_graph, FT, FL, Fg, FT+FL+Fg);

		  //q2charm[i], xmod, sigma(xmod, q2charm[i],0.0,par_temp));
		  //q2charm[i], x_Bj, sigma(xmod, q2charm[i],0.0,par_temp));
		  //q2charm[i], x_Bj, xmod);
	}
	fprintf(fout,"\n"); 
    }

    fclose(fout);

    //update_graph("f_2_charm.plt");
}

/*******************************************************************************
* The function plotting 
*******************************************************************************/
void graph_diff_zeus94 (void) {

    FILE *fout; 
     int i,k;

    double q2diff[4] = {8.0,14.0,27.0,60.0};
    double dbeta_graph = 0.0, FT = 0.0,FL = 0.0 ,Fg;
    //double par_temp[5];

    double npoints; /* Number of points used to plot */

    if (model == 1)
	chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);

    xp = 0.0042;
    xmod = xp;

    npoints = 10;

    fout = fopen("graph_diff_zeus94.dat","w");    

    for(i=0;i<=3;i++) {
	for(k=0;k<npoints;k++) {
	    dbeta_graph = 0.105*k+0.01; /* x = 10^-1...10^-6 */

	    FT = FD_T(q2diff[i],dbeta_graph,xp,m_fsq);
	    FL = FD_L(q2diff[i],dbeta_graph,xp,m_fsq);
	    Fg = FD_gluon(q2diff[i],dbeta_graph,xp);

	    fprintf(fout,"%e  %e  %e  %e  %e\n", 
			  dbeta_graph, FT, FL, Fg, FT+FL+Fg);
	}
	fprintf(fout,"\n"); 
    }
    fclose(fout);
}

/*******************************************************************************
* The function plotting 
*******************************************************************************/
void graph_diff_xp(void) {

    FILE *fout; 
     int i,k,j;

    double q2diff[7] = {2.7,4,6,8,14,27,55};
    double dbeta_graph[7][6]= {{0.652,0.231,0.0698,0.0218,0.0067,0.0030},\
                               {0.735,0.308,0.1000,0.0320,0.0099,0.0044},\
      			       {0.807,0.400,0.1430,0.0472,0.0148,0.0066},\
			       {0.848,0.471,0.1820,0.0620,0.0196,0.0088},\
			       {0.907,0.609,0.2800,0.1040,0.0338,0.0153},\
			       {0.949,0.750,0.4290,0.1820,0.0632,0.0291},\
			       {0.975,0.859,0.6040,0.3130,0.1210,0.0000}};

    double dxp,FT,FL,Fg = 0.0;

                         /*{0.949,0.750,0.4290,0.1820,0.0632,0.0291},\ */
    double npoints; /* Number of points used to plot */

    if (flavour == 1)
	chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);


    npoints = 4;

    fout = fopen("dxp_model.dat","w");    
    //fout = fopen("diff_xpMODEL.dat","w");    

    for(k=0;k<6;k++) {     /* Column, 6 */
	for(i=0;i<7;i++) {         /* Row , 7 */
	    for(j=0;j<npoints;j++) {
		//dxp = 0.01;
		dxp = 0.05*pow(10,-j);
		xmod = dxp;
		FT = FD_T(q2diff[i],dbeta_graph[i][k],xp,m_fsq);
		FL = FD_L(q2diff[i],dbeta_graph[i][k],xp,m_fsq);
		Fg = FD_gluon(q2diff[i],dbeta_graph[i][k],xp);

		//printf("%e  %f\n",dxp, dbeta_graph[i][k]);
		printf("%4.1f  %4.2f     %e  %e  %e\n",
                        q2diff[i], dbeta_graph[i][k],dxp, FT+FL,Fg);
		fprintf(fout,"%e  %e  %e\n",dxp, FT+FL,Fg);
            }
	    fprintf(fout,"\n"); 
	}
	//printf("Row %d finished \n",i+1); 
	printf("Column %d finished\n",k+1); 
    }
    fclose(fout);
}

/*******************************************************************************
* The function plotting 
*******************************************************************************/
void graph_diff_xp_9(void) {

    FILE *fout; 
     int i,k,j;

    double q2diff[3] = {2.7,8.0,27.0};
    double dbeta_graph[3][3]= {{0.652,0.231,0.0218},\
                               {0.848,0.471,0.0620},\
                               {0.949,0.750,0.182}}; 
    double dxp,FT,FL,Fg = 0.0;

    double npoints; /* Number of points used to plot */

    if (flavour == 1)
	chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);


    npoints = 4;

    fout = fopen("diff_xpMODEL.dat","w");    

    for(i=0;i<2;i++) {         /* Row */
	for(k=0;k<3;k++) {     /* Column */
	    for(j=0;j<npoints;j++) {
		//dxp = 0.01;
		dxp = 0.05*pow(10,-j);
		xmod = dxp;
		FT = FD_T(q2diff[i],dbeta_graph[i][k],xp,m_fsq);
		FL = FD_L(q2diff[i],dbeta_graph[i][k],xp,m_fsq);
		Fg = FD_gluon(q2diff[i],dbeta_graph[i][k],xp);

		//printf("%e  %f\n",dxp, dbeta_graph[i][k]);
		printf("%e  %e\n",dxp, FT+FL+Fg);
		fprintf(fout,"%e  %e\n",dxp, FT+FL+Fg);
            }
	    fprintf(fout,"\n"); 
	}
    }
    fclose(fout);
}
