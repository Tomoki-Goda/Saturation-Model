/*  
 *      Name: gbw
 *    Author: Sebastian Sapeta
 *      Date: August 2005
 *
 *            This program is a general framework to perform calculations
 *            within GBW and BGK model. 
 */

#include     <stdio.h>
#include   "control.h"
#include "test.h"
#include  "cfortran.h"
#include    "minuit.h"
//#include    "diffraction.h"



int main(int argc, char *argv[]) {


double  arglist[10];

/******************************************************************************/
/* Fit GBS                                                                    */
/******************************************************************************/

if ((action == 0)&&(model==2)) { 

    int error_flag = 0;
    
    /* Print starnting information */
    printinfo("stdout");

    /* Read data to global arrays */ 
    if (datatype == 1) {
      readdata_new();
    } else {
      readdata();
    }
    /* Initialize Minuit */
    MNINIT(5,6,7);

    /* Parameters definition */
    MNPARM(1,"sigma_0",
	   par_start[0],par_error[0],par_min[0],par_max[0],error_flag);
    MNPARM(2," lambda",
	   par_start[1],par_error[1],par_min[1],par_max[1],error_flag);
    MNPARM(3,"    x_0",
	   par_start[2],par_error[2],par_min[2],par_max[2],error_flag);
    if(sudflag>=1){
        MNPARM(4,"     C",
	     par_start[3],par_error[3],par_min[3],par_max[3],error_flag);
        MNPARM(5,"   mu02",
	       par_start[4],par_error[4],par_min[4],par_max[4],error_flag);
    };
    if(sudflag==2){
	MNPARM(4,"    g1",
			par_start[5],par_error[5],par_min[5],par_max[5],error_flag);
       // MNPARM(7,"    g2",
         //      par_start[6],par_error[6],par_min[6],par_max[6],error_flag);
	
    }; 

    arglist[0] = STRategy;

    /* Set strategy to STRategy from main.h, 0-fast, 1-default, 2-precise */
    MNEXCM(fcn,"SET STR",arglist, 1,error_flag,0);

    /* Minimalization procedure MIGRAD */
    MNEXCM(fcn,"MIGRAD",0,0,error_flag,0);
    
    /*FILE * resfile=fopen("resultGBW.txt","w");
    fprintf(resfile, "sigma0 = %f \n lambda=%f \n x_0 = %f \n C = %f \n mu02 = %f \n ",sigma0,lambda,x_0,C,mu02 );
    fprintf(resfile, "Q_max = %f/n Q_min = %f/n x_max = %f/n ",q_up,q_down,x_up);
    fclose();*/
    
    
}



/******************************************************************************/
/* Fit GBW                                                                    */
/******************************************************************************/

if ((action == 0)&&(model==0)) { 

    int error_flag = 0;
    
    /* Print starnting information */
    printinfo("stdout");

    /* Read data to global arrays */ 
    if (datatype == 1) {
      readdata_new();
    } else {
      readdata();
    }
    /* Initialize Minuit */
    MNINIT(5,6,7);

    /* Parameters definition */
    MNPARM(1,"sigma_0",
	   par_start[0],par_error[0],par_min[0],par_max[0],error_flag);
    MNPARM(2," lambda",
	   par_start[1],par_error[1],par_min[1],par_max[1],error_flag);
    MNPARM(3,"    x_0",
	   par_start[2],par_error[2],par_min[2],par_max[2],error_flag);


    arglist[0] = STRategy;

    /* Set strategy to STRategy from main.h, 0-fast, 1-default, 2-precise */
    MNEXCM(fcn,"SET STR",arglist, 1,error_flag,0);

    /* Minimalization procedure MIGRAD */
    MNEXCM(fcn,"MIGRAD",0,0,error_flag,0);
}


/******************************************************************************/
/* Fit BGK                                                                    */
/******************************************************************************/

if ((action == 0)&&(model==1)) { 
   
    int     error_flag = 0;

    char *fitres_file = "fitres.res";
        /* Print starnting information */
    system("clear");   
    printinfo("stdout");
    printinfo(fitres_file);

    /* Read data to global arrays */ 
    if (datatype == 1) {
      readdata_new();
    } else {
      readdata();
    }

    /* Initialize Minuit */
    MNINIT(5,6,7);

    /* Parameters definition */
    MNPARM(1,"sigma_0 ",
	  bgk_parst[0],bgk_parer[0],bgk_parmin[0],bgk_parmax[0],error_flag);
    MNPARM(2,"A_g     ",
	  bgk_parst[1],bgk_parer[1],bgk_parmin[1],bgk_parmax[1],error_flag);
    MNPARM(3,"lambda_g",
	  bgk_parst[2],bgk_parer[2],bgk_parmin[2],bgk_parmax[2],error_flag);
    MNPARM(4,"C       ",
	  bgk_parst[3],bgk_parer[3],bgk_parmin[3],bgk_parmax[3],error_flag);
    MNPARM(5,"mu02    ",
	  bgk_parst[4],bgk_parer[4],bgk_parmin[4],bgk_parmax[4],error_flag);


    arglist[0] = STRategy;
        
    /* Set strategy to STRategy from main.h, 0-fast, 1-default, 2-precise */
    MNEXCM(fcn,"SET STR",arglist, 1,error_flag,0);

    arglist[0] = 500;
    arglist[1] = 0.1;

    /* Minimalization procedure MIGRAD */
    MNEXCM(fcn,"MIGRAD",arglist,2,error_flag,0);

    //fclose(res_file);
/*
    int num, ivarbl;
    char *par_name;
    double val[1], error[1], *bnd1, *bnd2;

    num = 1;
    MNPOUT(num,par_name,val[0],error[0],bnd1,bnd2,ivarbl);
    printf("sigma_0   %e  +/-  %e\n",val[0], error[0]);
    num = 2;
    MNPOUT(num,par_name,val[0],error[0],bnd1,bnd2,ivarbl);
    printf("A_g       %e  +/-  %e\n",val[0], error[0]);
    num = 3;
    MNPOUT(num,par_name,val[0],error[0],bnd1,bnd2,ivarbl);
    printf("lambda_g  %e  +/-  %e\n",val[0], error[0]);
    num = 4;
    MNPOUT(num,par_name,val[0],error[0],bnd1,bnd2,ivarbl);
    printf("C         %e  +/-  %e\n",val[0], error[0]);
    num = 5;
    MNPOUT(num,par_name,val[0],error[0],bnd1,bnd2,ivarbl);
    printf("mu02      %e  +/-  %e\n",val[0], error[0]);
*/
    graph_dipolcs();
    graph_f_2_charm();
    xgpdf_av();
}



/******************************************************************************/
/* List function values for all experimental data                             */
/******************************************************************************/

if (action == 1) { 

    int licznik;
    //double sxmod;
    //char *result_file = "wynik";
    //FILE *res_file;
    FILE *list_file;

    //res_file = freopen("res_rile.txt","w",stdout);
    list_file = fopen("f2_list.dat","w");

    /* Print starnting information */
    system("clear");   
    printinfo("stdout");
    //printinfo(result_file);

    /* Read data to global arrays */ 
    readdata();
    switch (model) {
    case 0:
	sigma_0  = par_start[0]; 
	lambda   = par_start[2];
	x_0      = par_start[1];
    break;
    case 1:
	/* Set parameters equal to start fit parameters */
	sigma_0  = bgk_parst[0]; 
	A_g      = bgk_parst[1]; 
	lambda_g = bgk_parst[2]; 
	C        = bgk_parst[3]; 
	mu02     = bgk_parst[4]; 
	/* Perform chebyshev approximation for xgpdf */
	chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);
    break;
    }

/*
    double spar[5]  = {15.233,3.2748,0.000,1.2083,1.1018};
    //double spar[5]  = {15.74,3.251,0.000,1.398,0.9055};
    double gpar[5]  = {15.65,3.908,0.000,2.313,0.8628};

    sigma_0  = spar[0]; 
    A_g      = spar[1]; 
    lambda_g = spar[2]; 
    C        = spar[3]; 
    mu02     = spar[4]; 

*/
    /*
    sigma_0  = gpar[0]; 
    A_g      = gpar[1]; 
    lambda_g = gpar[2]; 
    C        = gpar[3]; 
    mu02     = gpar[4]; 
    */

    //double sfit = 0.0, gfit = 0.0;
    //double schi = 0.0, gchi = 0.0;


    
    double chip_sum = 0.0;

    printf("%.10f\n",1.0/(pi*pi*pi*0.389379));
    printf("%.10f\n",1.0/(PI*PI*PI*0.389379));
    //for(licznik=284;licznik<NDATA;licznik++) {
    for(licznik=0;licznik<NDATA;licznik++) {
        
	xmod = xmoddata[licznik];
	Q2   = q2data[licznik];
	y    = ydata[licznik];

        /*
	printf("%e  %e  %e\n", xmod, Q2,sigma_bgk_cheb(1/sqrt(Q2)));
        */

	//sigmafcn[licznik] = sigma(xmod,Q2,y,par_start);
        /*
        temp_var = 0;
	sfcn0[licznik] = sigma(xmod,Q2,y,bgk_parst);
        temp_var = 1;
	sfcn1[licznik] = sigma(xmod,Q2,y,bgk_parst);

	printf("%e  %e  %e  %1.10e  %1.10e\n", 
               xmod, Q2,csdata[licznik],sfcn0[licznik],sfcn1[licznik]);
        */                           
        /*
        //double mu2 = 0.000001;
        double mu2 = mu02;
        //double xgpdf_cheb = sigma_bgk(0.001);
        double xgpdf_cheb = sigma_bgk(10000);
        */
/*
        double mu2 = Q2;
        double xgpdf_cheb = chebev(xmin,xmax,Qmin,Qmax,MX,MQ,coef,xmod,mu2); 
        printf("%.1e  %e %.10e \n",Q2,xmoddata[licznik],xgpdf_cheb);
*/
    /*    
	printf("%e  %e  %1.10e  %1.10e\n", 
               xmod, Q2, sigma_bgk(Q2),fabs(sigma_bgk(Q2)-sigma_bgk_cheb(Q2)));
    */
	
       /* 

        double sxb,xper,sper;
        sxb   = sigma(xbdata[licznik],Q2,y,bgk_parst);
        sxmod = sigma(xmoddata[licznik],Q2,y,bgk_parst);
        xper  = xmoddata[licznik]/xbdata[licznik]*100;
        sper  = sxmod/sxb*100;
        double a[2] = {1.0e-17,1.0e-17};
        double b[2] = {Rmax,0.5};

        printf("                                                   %.10e\n",
               simps2d(a,b,1.0e-06,&uif_bgk));
	printf("%.1e  %e  %e  %e  %.10e  %.2f  %.2f\n",Q2,
               //xbdata[licznik],xmoddata[licznik],sxb,sxmod,xper,sper);
               xbdata[licznik],xmoddata[licznik],sxb,sxmod);
        */


        /*
        schi = 
        (csdata[licznik]-sigma(xmoddata[licznik],Q2,y,spar))/erdata[licznik];

        gchi = 
        (csdata[licznik]-sigma(xmoddata[licznik],Q2,y,gpar))/erdata[licznik];

        sfit +=schi*schi;
        gfit +=gchi*gchi;

        */

/*
	printf("%e  %e  %1.10e  %1.10e\n", 
              //xmod, Q2, sigma_bgk(Q2),fabs(sigma_bgk(Q2)-sigma_bgk_cheb(Q2)));
               //xmod, Q2, sigma_bgk(Q2),sigma_bgk_cheb(Q2));
               xmod, Q2, xgpdf(xmod,Q2),sigma_bgk_cheb(Q2));
*/
/*
        printf("%3d  %.1e  %e  %e %.10e   %.10e\n",
               licznik+1, Q2, xbdata[licznik],xmoddata[licznik],
               csdata[licznik],
               //sigma(0.02,Q2,y,spar));
               sigma(xmoddata[licznik],Q2,y,bgk_parst));
*/
    double chip;
    switch (model) {
    case 0:
	sigmafcn[licznik]=sigma(xmoddata[licznik],
			q2data[licznik],ydata[licznik],par_start);
    break;
    case 1:
	sigmafcn[licznik]=sigma(xmoddata[licznik],
			q2data[licznik],ydata[licznik],bgk_parst);
    break;
    }


    chip = (csdata[licznik]-sigmafcn[licznik])/erdata[licznik];
    chip_sum += chip*chip;

    printf("%3d  %6.2f  %f  %f  %f  %f  %8.2f %8.2f\n",
                licznik+1,q2data[licznik], xmoddata[licznik],erdata[licznik],
                csdata[licznik],sigmafcn[licznik],chip*chip,chip_sum);

    fprintf(list_file,"%3d  %6.2f  %f  %f  %f  %f  %8.2f %8.2f\n",
                licznik+1,q2data[licznik], xmoddata[licznik],erdata[licznik],
                csdata[licznik],sigmafcn[licznik],chip*chip,chip_sum);
    //if(((licznik+1)==58)||((licznik+1)==102)||((licznik+1)==172))
      //	fprintf(list_file,"\n");
    }
    
    fclose(list_file);
    //fclose(res_file);
    //graph_dipolcs();
    //graph_f_2_charm();
    //xgpdf_av();
    //graphdata();
}


/******************************************************************************/
/* Testing module              					              *//******************************************************************************/
if (action == 2) { 

    int test_type = 0;

    system("clear");   /* Clear sreen */
    printf("-----------------------------------------------------------------");
    printf("--------------\n");
    printf("                               TESTING MODULE\n");
    printf("-----------------------------------------------------------------");
    printf("--------------\n\n");
    printf("        1 - compare my xpdf with K.Golec-Biernat xpdf\n\n");
    printf("        2 - compare exact xpdf with xpdf_cheb        \n\n");
    printf("        3 - compare red_cs calculated with dadmul and simps2d\n\n");
    printf("        4 - test dadmul integration function\n\n");
    printf("        5 - psi charm test\n\n");
    printf("        ");
    scanf("%d",&test_type);

    switch (test_type) {
        case 1:
	comp_xgpdf_mg();
        break;
        case 2:
	comp_xgpdf_ec();
        break;
        case 3:
	comp_int();
        break;
        case 4:
	dadmul_test();
        break;
        case 5:
	psi_charm_test();
        break;

    }
}



/******************************************************************************/
/* Something else, not yet finished                                          */
/******************************************************************************/
if (action == 100) { 

    int ln;       

    /* Read data to global arrays */ 
    readdata();

    sigma_0  = bgk_parst[0]; 
    A_g      = bgk_parst[1]; 
    lambda_g = bgk_parst[2]; 
    C        = bgk_parst[3]; 
    mu02     = bgk_parst[4]; 


    chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);

    xmod = 1.0e-02;

    for(ln=0;ln<16;ln++)
    printf("BWK_C:  %e  %e \n", xmoddata[ln], sigma_bgk_cheb(ln*0.5+1.0));
    //printf("BWK:    %e  %e \n", xB, sigma_bgk(10.0));


    sigma_0 = 23.0; 
    lambda  = 0.29; 
    x_0     = 3.0e-04; 

    //printf("GBW:  %e  %e \n", xB, sigma_gbw(1.0));


    //printf("%e  %e \n", xB, xgpdf(xB,10));
    /* List function values */
/*
    for(i=0;i<0;i++) {
	xB = xmoddata[i];
	printf("%e  %e \n", xB, sigma_bgk(0.01*i+0.5));
    }
*/
}
if (action == 4) { 
/*
    sigma_0  = bgk_parst[0]; 
    A_g      = bgk_parst[1]; 
    lambda_g = bgk_parst[2]; 
    C        = bgk_parst[3]; 
    mu02     = bgk_parst[4]; 
*/

    sigma_0 = par_start[0]; 
    lambda  = par_start[1]; 
    x_0     = par_start[2]; 

 //   chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);

//double bessk(int n,double x); 

    dbeta = 0.209;
    Q2 = 14.0;
    xp = 0.0042;
    xmod = xp;
    //double ttt = 0.8;
    //double ztt = 0.2;
    //double fd;
    //printf("% f\n",bessk(2,1.2) );
    //printf("% f\n",bessj(2,1.2) );
    //printf("% f\n",phi_2_integrand(&ttt));
    //printf("% f\n",phi_2(&ttt));
    //printf("% f\n",FD_gluon_integrand(&ztt));
    //printf("% f\n",FD_gluon(14, 0.2,0.0042 ));
    //printf("% f\n",FD_gluon(14, 0.3,0.0042 ));
    //printf("% f\n",FD_gluon(14, 0.4,0.0042 ));
    //printf("% f\n",FD_gluon(Q2, dbeta, xp));
    //printf("% f\n",phi_2(&ttt));
    xmod = xp;
    //printf("% f\n",FD_gluon(14, 0.2,xp ));
    xmod = xp;
    //printf("% f\n",FD_T(8, 0.2, xp, m_fsq));

/*
    printf("% f\n",phi_0(&ttt) );
    printf("% f\n",phi_1(&ttt) );
    printf("% f\n",FD_L_integrand(&ttt));
    printf("% f\n",FD_L(Q2, dbeta, xp, m_fsq));
    printf("% f\n",FD_T(Q2, dbeta, xp, m_fsq));

    fd =   FD_T(Q2, dbeta, xp, m_fsq) + FD_T(Q2, dbeta, xp, m_ch)
         + FD_L(Q2, dbeta, xp, m_fsq) + FD_L(Q2, dbeta, xp, m_ch);

    printf("suma:\n% f \n",fd);
*/
}

if (action == 5) { 

    model = 1;
    switch (model) {
    	case 0:
	sigma_0 = par_start[0]; 
	lambda  = par_start[1]; 
	x_0     = par_start[2]; 
	break;
        case 1:
	sigma_0  = bgk_parst[0]; 
	A_g      = bgk_parst[1]; 
	lambda_g = bgk_parst[2]; 
	C        = bgk_parst[3]; 
	mu02     = bgk_parst[4]; 
	break;
    }

double alpha_temp(double L, int nf) {

    double b0;
    b0 = (33-2*nf)/(12*pi);
    return 1/(b0*log(91*91/(L*L*1.0e-06))); 
}

/*
printf("l_qcd  nf   alpha_s\n");
printf("%.0f    %d    %.4f\n",200.0,0, alpha_temp(200.0,0));
printf("%.0f    %d    %.4f\n",250.0,0, alpha_temp(250.0,0));
printf("%.0f    %d    %.4f\n",300.0,0, alpha_temp(300.0,0));
printf("%.0f    %d    %.4f\n",200.0,5, alpha_temp(200.0,5));
printf("%.0f    %d    %.4f\n",250.0,5, alpha_temp(250.0,5));
printf("%.0f    %d    %.4f\n",300.0,5, alpha_temp(300.0,5));

chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);
double q2_photo;
double norm_photo;

double sigma_photo;
norm_photo = 4*PI*PI/(137*2.568192);
q2_photo = 1.3e+00;
*/

//sigma_photo = norm_photo*1000*sigma_c(0.0001,q2_photo,0,bgk_parst)/q2_photo;
//printf("%f\n",sigma_photo);



/*

q2_photo = 0.0;
sigma_photo = norm_photo*sigma_c(0.01,q2_photo,0,bgk_parst);
printf("%f\n",sigma_photo);
sigma_photo = norm_photo*sigma_b(0.01,q2_photo,0,bgk_parst);
printf("%f\n",sigma_photo);


*/


/*
q2_photo = 1.0e-03;
sigma_photo = norm_photo*sigma_c(0.01,q2_photo,0,bgk_parst)/q2_photo;
printf("%f\n",sigma_photo);
q2_photo = 1.0e-04;
sigma_photo = norm_photo*sigma_c(0.01,q2_photo,0,bgk_parst)/q2_photo;
printf("%f\n",sigma_photo);
q2_photo = 1.0e-05;
sigma_photo = norm_photo*sigma_c(0.01,q2_photo,0,bgk_parst)/q2_photo;
printf("%f\n",sigma_photo);
q2_photo = 1.0e-06;
sigma_photo = norm_photo*sigma_c(0.01,q2_photo,0,bgk_parst)/q2_photo;
printf("%f\n",sigma_photo);
q2_photo = 1.0e-07;
sigma_photo = norm_photo*sigma_c(0.01,q2_photo,0,bgk_parst)/q2_photo;
printf("%f\n",sigma_photo);
*/

    printf("%.8f\n",1.0/(pi*pi*pi));
    //readdata();
    //count_data();
    //graph_dipolcs();
    //graph_dipole_scaling();
    //graph_sat_scale();
    //graph_sat_scale_05();
    //graph_f_2_beauty(0);
    //graph_f_2_beauty(1);
    //graph_flq2(0);
    alpha_fpdf_graph();
    //graph_flq2(0,0);
    //graph_flq2(0,1);
    //graph_flq2(1,1);
    //graph_fl();
    //graph_alphagluon();
    //graph_f_2_charm_9(0);

    //graph_f2c9q2(0);
    //graph_f2c9q2(1);
    //graph_f2b9q2(0);
    //graph_f2b9q2(1);

    //graph_f2chf2();
    //graph_gluon();

    //graph_f_2_charm_9(1);

    //graph_f2ch_comp();

    //graph_diff_xp();
    //graph_diff_xph1zeus();
    //graph_diff();
    //graph_diff_zeus94();
    //printf("%f\n", xgpdf_av());
    //w2_test();
    //xmod = 0.0001;
    //printf("%f \n",FD_T (10, 0.5, 0.0001, 0.0));

    //printf("%f\n",sigma_DT(100.0, 100.0, 3.0, 4.0, 0.0));
    //printf("%f\n", find_zero (0,5, &prob));
    //graph_crit_line(1);

/*
    chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);
    xmod = 1.0e-06;
    printf("%f\n", find_zero (0,50, &crit_func));
*/

/*
    chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);
    double _xmod = 0.000554;
    double _Q2 = 35;
   printf("%f\n",sigma_l(_xmod,_Q2,y,bgk_parst)+sigma_c(_xmod,_Q2,y,bgk_parst));
*/
}
return 0;
}

