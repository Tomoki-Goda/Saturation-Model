/*  
 *
 */

#include  <math.h>
#include <string.h>
#include "gluons.h"
#include "float.h"
#include "chebyshev.h"
#include "chebyshev3.h"
//#include"photon-wave-function.h"
//#include "cuba.h"
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_spline.h>

#define NTRY 100

double arg_unint_x;
double arg_unint_l2;
double rp;
double x4int;
double q24int;
double R4int;


extern double dzerox_();

/*******************************************************************************
*******************************************************************************/
/*double determine_rmax(double (*func)(double)) {
  int imax = 20;
  for (int i = 0; i<imax; i++) {
    double rcut = rmin + (rmax-rmin)*i/imax;
    printf("RMAX: %e %e\n",rcut, func(rcut));
  }
}
*/
/*******************************************************************************
* The 'light' only photon wave function - reduction cross section and F_2 form
*******************************************************************************/
double psisql (double r, double z)  {   

    double     z_bar =  z*z+(1-z)*(1-z);
    double     y_bar =  (y*y)/(1+(1-y)*(1-y));
    double   Qsq_bar =  z*(1-z)*Q2+m_fsq;
    double      Qsq2 =  sqrt(Qsq_bar)*r;
    double        nl;

    double bessel_k0 = dbesk0_(&Qsq2);
    double bessel_k1 = dbesk1_(&Qsq2);

    double value;
   
    //nl = 1.0;
    nl = 5.0/6.0;
    switch(dataform) {
    case 0: /* Reduced cross section form */  
	value =  norm * nl * Q2*(z_bar*Qsq_bar*bessel_k1*bessel_k1
                             + (m_fsq+(1-y_bar)*4*Q2*z*z*(1-z)*(1-z))
                                  *bessel_k0*bessel_k0);
        break;
    case 1: /* F_2 form */  
	value =  norm * nl * Q2*(z_bar*Qsq_bar*bessel_k1*bessel_k1
                                 + (m_fsq+4*Q2*z*z*(1-z)*(1-z))
                                 *bessel_k0*bessel_k0);
 	//value =  norm * nl * Q2;	
        break;
    case 2: /* F_L form */
	value =  norm * nl * Q2*4*Q2*z*z*(1-z)*(1-z)*bessel_k0*bessel_k0;
        break;
    }

   return value;
}

/*******************************************************************************
* The 'strange' only photon wave function - reduced cross section and F_2 form
*******************************************************************************/
double psisqs (double r, double z)  {   

    double     z_bar =  z*z+(1-z)*(1-z);
    //double     y_bar =  (y*y)/(1+(1-y)*(1-y));
    double   Qsq_bar =  z*(1-z)*Q2+m_s;
    double      Qsq2 =  sqrt(Qsq_bar)*r;
    double       ns;

    double bessel_k0 = dbesk0_(&Qsq2);
    double bessel_k1 = dbesk1_(&Qsq2);

    double value;
   
/*
	value =  norm * Q2*4*Q2*z*z*(1-z)*(1-z)*bessel_k0*bessel_k0;
*/
    //ns = 1.0;
    ns = 1.0/6.0;

    switch (dataform) {
    case 1: /* F_2 form */  
	value = norm*ns*Q2*(z_bar*Qsq_bar*bessel_k1*bessel_k1
                + (m_s+4*Q2*z*z*(1-z)*(1-z))*bessel_k0*bessel_k0);
	break;
    case 2: /* F_L form */
	value = norm*ns*Q2*4*Q2*z*z*(1-z)*(1-z)*bessel_k0*bessel_k0;
        break;
    }

   return value;
}
/*******************************************************************************
* The 'charm' only photon wave function - reduction cross section and F_2 form
*******************************************************************************/
double psisqc (double r, double z)  {   

    double     z_bar =  z*z+(1-z)*(1-z);
    //double     y_bar =  (y*y)/(1+(1-y)*(1-y));
    double   Qsq_bar =  z*(1-z)*Q2+m_ch;
    double      Qsq2 =  sqrt(Qsq_bar)*r;
    double norm_c;

    double bessel_k0 = dbesk0_(&Qsq2);
    double bessel_k1 = dbesk1_(&Qsq2);

    double value;
   
/*
	value =  norm * Q2*4*Q2*z*z*(1-z)*(1-z)*bessel_k0*bessel_k0;
*/
    switch (dataform) {
    case 1: /* F_2 form */  
	value = (2.0/3.0)*norm*Q2*(z_bar*Qsq_bar*bessel_k1*bessel_k1
                + (m_ch+4*Q2*z*z*(1-z)*(1-z))*bessel_k0*bessel_k0);
 	//value = (2.0/3.0)*norm * Q2;	
	break;
    case 2: /* F_L form */
	value = (2.0/3.0)*norm * Q2*4*Q2*z*z*(1-z)*(1-z)*bessel_k0*bessel_k0;
        break;
    }

    if (photo ==1) {
        norm_c = (2.0/3.0)*(2.0/3.0);
	value = norm_c*z_bar*Qsq_bar*bessel_k1*bessel_k1;
    }


   return value;
}

/*******************************************************************************
* The 'bottom' only photon wave function - reduction cross section and F_2 form
*******************************************************************************/
double psisqb (double r, double z)  {   

    double     z_bar =  z*z+(1-z)*(1-z);
    double   Qsq_bar =  z*(1-z)*Q2+m_b;
    double      Qsq2 =  sqrt(Qsq_bar)*r;

    double bessel_k0 = dbesk0_(&Qsq2);
    double bessel_k1 = dbesk1_(&Qsq2);

    double value;

    switch (dataform) {
    case 1: /* F_2 form */  
	value = (1.0/6.0)*norm*Q2*(z_bar*Qsq_bar*bessel_k1*bessel_k1
                + (m_b+4*Q2*z*z*(1-z)*(1-z))*bessel_k0*bessel_k0);
 	//value = (1.0/6.0)*norm * Q2;	
	break;
    case 2: /* F_L form */
	value = (1.0/6.0)*norm * Q2*4*Q2*z*z*(1-z)*(1-z)*bessel_k0*bessel_k0;
        break;
    }

    if (photo ==1) {
	value = (1.0/6.0)*norm*z_bar*Qsq_bar*bessel_k1*bessel_k1;
    }


   return value;
}
/*******************************************************************************
* The 'light+charm' photon wave function - reduction cross section and F_2 form
*******************************************************************************/
double psisqlc (double r, double z)  {   

    double value;
    value =  psisql(r,z)+psisqc(r,z);  
    return value;
}

/*******************************************************************************
* The dipole-proton GBW cross section 
*******************************************************************************/
double sigma_gbw (double r) {

   return sigma_0*(1-exp(-0.25*r*r*pow(x_0/xmod,lambda)));
}

/*******************************************************************************
* 
*******************************************************************************/
double sudakov(double r, double mu2) {
    
    //if (1/r*r < Lambda2) return 0.0; 
    //double bmax2 = (1.26095)/(C*C);
    //double bmax2 = pow(bmax,2.0);

    //double mub2 = 1.26095/(r*r/(1+r*r/bmax2));
    //double mub2=C/(pow(r,2)) + mu02;
    double mub2=C*(1/pow(r,2)+1/pow(r_max,2));
    //double mub2=C*(1/(pow(r,2)) + 1/pow(bmax,2));
    //printf("sud %e %e %e\n", r, mub2, mu2);
    //double mub2 = 1.26095/(r*r);
    if (mu2 < Lambda2 || mub2 < Lambda2|| mu2 < mub2) {/*printf("1"); */return(0.0);}; 
    // if (pow(r,2.0) > bmax2){/*printf("2");*/ return(0.0);}; 
       

    //printf("%e %e %e\n", r, mub2, mu2);
    double CA = 3.0;
    double b0 = (33-2*n_f)/12;
   // double b0 = (33-2*n_f)/(12*pi);
    
    //printf("%e\n", sqrt(bmax2)/sqrt(bmax2*mu2-1));

    //if (r>sqrt(bmax2)/sqrt(bmax2*mu2-1)) return 0.0;
    //printf("sud");
    double val = CA/(2*b0*PI)*
                 (log(mu2/Lambda2)*log(log(mu2/Lambda2)/log(mub2/Lambda2)) - 
	         log(mu2/mub2));
    //double nll = -CA/PI*log(log(mu2/Lambda2)/log(mub2/Lambda2));
    //val += nll;
    //val  = val/exp(mu02*r);
    return val;
}

/******************************************************************************
 *non perturbative sudakov...
 * ***************************************************************************/
double sudakov_np(double  r,double mu2){

	double val=g1 * pow(r,2.0)/(2.0)+ g2 * ( log(mu2/pow(Q0,2.0)) * log(1+pow(r/r_max,2)/*mu02*pow(r,2 )/C*/)/4.0 );
	//double val=g1 * pow(r,2.0)/(2.0) + g2 * ( log(mu2/pow(Q0,2.0)) * log( 1+(mu02* pow(r,2)/C ))/4.0 );
	//double val=g1 * pow(r,2.0)/(2.0) + g2 * ( log(mu2/pow(Q0,2.0)) * log( 1+(pow(r,2)/pow(bmax,2) ) )/4.0 );

	//double val=g1 * pow(r,2.0)/(2.0) + g2 * ( log(mu2/pow(Q0,2.0)) * log((pow(bmax,2.0)+pow(r,2.0))/pow(bmax,2.0))/4.0 );

	return(val);
}

/*******************************************************************************
* 
*******************************************************************************/
double sudakov2(double r, double mu2) {

    //double bmax2 = C*C;
    //double mub2 = 1.26095/(r*r/(1+r*r/bmax2));
    double mub2=C/(pow(r,2)) + mu02;
    if (mu2 < mub2 || mu2 < Lambda2 || mub2 < Lambda2) return 0.0; 

    double CA = 3.0;
    double b0 = (33-2*n_f)/(12*pi);
    double as = 0.2;
    double val = as*CA*pow(log(Q2/mub2),2)/(4*pi);
    //double nll = - as*b0*CA*log(Q2/mub2)/pi;
    //val += nll;

    return val;
}

/*******************************************************************************
* 
*******************************************************************************/
double S_gbs_int(double* r) {

  double val;
  double Qs2 = pow(x_0/x4int,lambda);
  //double gbw = r[0]*exp(-Qs2*r[0]*r[0]/4)*Qs2*(Qs2*r[0]*r[0]-4)/4;
  //double gbw = r[0]*log(r[0])*exp(-Qs2*r[0]*r[0]/4)*Qs2*(Qs2*r[0]*r[0]-4)/4;
  double gbw = 
         -r[0]*log(R4int/r[0])*exp(-Qs2*r[0]*r[0]/4)*Qs2*(Qs2*r[0]*r[0]-4)/4;

  if (sudflag >= 1) {
    double sud = exp(-sudakov(r[0],q24int));
    //double sud = exp(-sudakov(r[0],q24int));
    //printf("#%e %e %e\n", r[0], q24int, sud);
    
    val = gbw*sud;
    if(sudflag==2){
    val*=exp(-sudakov_np(r[0],q24int));
    }
  } else if (sudflag == 0) {
    val = gbw;
  }
  
  return val;
}


/*******************************************************************************
* x and q2 are passed to the integrant via global variables, as well as the
* lower limit of integration R
*******************************************************************************/
double S_gbs(double x, double q2, double r) {
  
  double rmi = 1.0e-06;
  double rma = 1.0e+03;
  double eps = 1.0e-06;

  extern double simps_();
  double dum1,dum2,dum3;
  double int_result;
  R4int = r;
  x4int= x;
  q24int = q2;
  //q24int=C/(pow(r,2.0))+mu02;

  //printf("%e %f %e\n", x, q2, r);
  simps_(&rmi,&r,&eps,&eps,&eps, &S_gbs_int,&dum1,&int_result,&dum2,&dum3);
  //simps_(&r,&rma,&eps,&eps,&eps, &S_gbs_int,&dum1,&int_result,&dum2,&dum3);

  //printf("%e %e %e\n", x, r, int_result);
  //printf("%e %e %e %e\n", x, q2, r, int_result);
  return int_result;
}

/*******************************************************************************
* 
*******************************************************************************/
double logS_gbs(double x, double q2, double r) {
  double val = S_gbs(x, q2, r);
  //if(val<DBL_MIN) val = log(DBL_MIN);
  //else val = log(val);
  //return val;
  return log(val);
}

/*******************************************************************************
* 
*******************************************************************************/
double S_gbs_cheb(double x, double q2, double r) {
  
   //printf("%e %e %e %e\n", xmod, q2, r,
   //       exp(chebev3(xmin,xmax,Qmin,Qmax,rmin,rmax,NX,NQ,NR,coef3,x,q2,r)));
   return exp(chebev3(xmin,xmax,Qmin,Qmax,rmin,rmax,NX,NQ,NR,coef3,x,q2,r));
}

/*******************************************************************************
* At this point, Q2 and xmod are set in sigma_l(). This function is not involved
* in the procedure of interpolation.
*******************************************************************************/
double sigma_gbs_exact(double r) {

   return sigma_0*S_gbs(xmod, Q2, r);
}

/*******************************************************************************
* At this point, Q2 and xmod are set in sigma_l(). This function is not involved
* in the procedure of interpolation.
*******************************************************************************/
double sigma_gbs (double r) {

   //printf("%e %e %e %e\n", xmod, Q2, r, S_gbs(xmod, Q2, r));
   //printf("%e %e %e %e %e\n", xmod, Q2, r, 
   //        S_gbs(xmod, Q2, r), S_gbs_cheb(xmod, Q2, r));
   return sigma_0*S_gbs(xmod, Q2, r);
   //return sigma_0*S_gbs_cheb(xmod, Q2, r);

   //if (r<0.5) return sigma_0*S_gbs_cheb(xmod, Q2, r);
   //else return sigma_gbw(r);

   //if (r<3) return sigma_0*S_gbs_cheb(xmod, Q2, r);
   //else sigma_0*S_gbs(xmod, Q2, r);
}

/*******************************************************************************
* The strong running coupling alpha_s
*******************************************************************************/
double alpha_s (double Q) {

    double b0;

    b0 = (33-2*n_f)/(12*pi);
    
    return 1/(b0*log(Q/Lambda2));  /* Q means Q^2 [GeV^2]*/

}
/*******************************************************************************
* The dipole-proton BGK cross section 
*******************************************************************************/
double sigma_bgk (double r) {
   
   //double mu2 = mu02/(1-exp(-mu02*r*r/C));
   double mu2 = C/(r*r) + mu02;
   double exnum; /* The numerator in the exp function */
   exnum = 0.389379*pi*pi*r*r*alpha_s(mu2)*xgpdf(xmod,mu2);

   //cacf = 9.0/4.0; 
   //return 8.0/3.0*sigma_0*(1-exp(-exnum*cacf/(3*sigma_0)));
   return sigma_0*(1-exp(-exnum/(3*sigma_0)));
}

/*******************************************************************************
* The dipole-proton BGK cross section  - Chebeshev approximation
*******************************************************************************/
double sigma_bgk_cheb (double r) {
   
   double mu2 = C/(r*r) + mu02;
   //double mu2 = mu02/(1-exp(-mu02*r*r/C));
   double exnum; /* The numerator in the exp function */
   double xgpdf_cheb;

   /*** This should be uncommented if sigma_x is used  ***/
   /*
   switch (light_charm) {
       case 0:
       xgpdf_cheb = chebev(xmin,xmax,Qmin,Qmax,MX,MQ,coef,xmod,mu2); 
       break;
       case 1:
       xgpdf_cheb = chebev(xmin,xmax,Qmin,Qmax,MX,MQ,coef,xmod1,mu2); 
       break;
  }
  */  
   xgpdf_cheb = chebev(xmin,xmax,Qmin,Qmax,MX,MQ,coef,xmod,mu2); 


   exnum = 0.389379*pi*pi*r*r*alpha_s(mu2)*xgpdf_cheb;
  
   /* If the xgpdf is negative print the warning message */
   if (exnum <0) {
   	negglu = 1;
        //printf("!!! %e   %e\n", xmod, mu2);
   }
   else
   	negglu = 0;

   //return sigma_0*(1.0-exp(-r*r));
   return sigma_0*(1.0-exp(-exnum/(3.0*sigma_0)));
}

/*******************************************************************************
* The dipole-proton BGK cross section 
*******************************************************************************/
double sigma_scaling (double r, double Q2_scal) {
   
   double exnum; /* The numerator in the exp function */
   exnum = r*r*Q2_scal;

   return sigma_0*(1-exp(-exnum));
}
/*******************************************************************************
* The critical function
*******************************************************************************/
double crit_func (double Q2s[0]) {
   
   double mu2 = C*Q2s[0]/4.0 + mu02;
   double xgpdf_cheb;
   double value; 

   xgpdf_cheb = chebev(xmin,xmax,Qmin,Qmax,MX,MQ,coef,xmod,mu2); 
   //value = 4*0.389379*pi*pi*alpha_s(mu2)*xgpdf_cheb/(3*sigma_0*Q2s[0])+log(0.3);
   value = 4*0.389379*pi*pi*alpha_s(mu2)*xgpdf_cheb/(3*sigma_0*Q2s[0])-1.0;
   //value = 4*0.389379*pi*pi*alpha_s(mu2)/(3*sigma_0*Q2s[0])-1.0;
   //p(-0.25*r*r*pow(x_0/xmod,lambda)));
  
   return value;
}

double crit_func_frozen (double Q2s[0]) {
   
   double mu2 =  mu02;
   double xgpdf_cheb;
   double value; 

   xgpdf_cheb = chebev(xmin,xmax,Qmin,Qmax,MX,MQ,coef,xmod,mu2); 
   value = 4*0.389379*pi*pi*alpha_s(mu2)*xgpdf_cheb/(3*sigma_0*Q2s[0])-1.0;
  
   return value;
}

double prob (double x[0]) {

    return x[0]*x[0]-4;
}

/*******************************************************************************
* Change x1 until f(x1)*f(x2) < 0 
*******************************************************************************/
int zbrac(double (*func)(double *), double *x1, double *x2) {

    int j;
    double f1,f2;

    if (*x1 == *x2) 
	printf("Bad initial range in zbrac \n");
    for (j=1;j<=NTRY;j++) {
    f1=(*func)(x1);
    f2=(*func)(x2);
	    if (f1*f2 < 0.0) {
                    return 1;
                    break;
            }
            else
		    *x1 += 0.01;
/*
	    else if (fabs(f1) < fabs(f2)) {
		    // *x1 += 0.1;
		    // *x1 += FACTOR*(*x1-*x2);
		    f1=(*func)(x1);
            }
	    else {
		    *x1 += 0.1;
		    // *x2 += FACTOR*(*x2-*x1);
                    //printf(":)\n");
		    f2=(*func)(x2);
            }
*/
    }
    return 0;
}

/*******************************************************************************
* Find zero of the function in the range (a,b) 
*******************************************************************************/
double find_zero (double a, double b, double (*func)(double *)) {

    double eps = 1.0e-10;
    int max = 1000;
    int mode = 2;
    double value;

    //printf("%e   %e\n",a,b);
    zbrac(func, &a, &b);
    //printf("%e   %e\n",a,b);
    value = dzerox_(&a,&b,&eps, &max, func, &mode);

    return value;

}
/*******************************************************************************
* The dipole-proton BGK cross section modified  - Chebeshev approximation
*******************************************************************************/
double sigma_bgk_cheb_mod (double r) {

   double mu2 = C/(r*r) + mu02;
   //double mu2 = mu02/(1-exp(-mu02*r*r/C));
   double exnum,cacf; /* The numerator in the exp function */
   double xgpdf_cheb;

   xgpdf_cheb = chebev(xmin,xmax,Qmin,Qmax,MX,MQ,coef,xmod,mu2); 
   exnum = 0.389379*pi*pi*r*r*alpha_s(mu2)*xgpdf_cheb;
  
   /* If the xgpdf is negative print the warning message */
   if (exnum <0) 
   	printf("!!! %e   %e\n", xmod, mu2);

   cacf = 9.0/4.0; 
   return sigma_0*(1-exp(-exnum*cacf/(3*sigma_0)));
}

/*******************************************************************************
* The integrand of DIS GBW cross section
*******************************************************************************/
extern double uif_gbs (int *dim, double *x) {

    /*  r =  x[0]    
        z =  x[1]  
     */

    double value;


    switch (flavour) {
    case 0: /* Light  only */
        value = x[0]*psisql(x[0],x[1])*sigma_gbs(x[0]);
        break;
    case 1: /* Light + charm */
        value = x[0]*psisqlc(x[0],x[1])*sigma_gbs(x[0]);
        break;
    case 2:
        break;
    default:
        printf("Error: wrong 'flavour' value");
        break;
    }

   return value;
}

extern double uif_gbs_l (int *dim, double *x) {

    /*  r =  x[0]    
        z =  x[1]  */

    double value;
    //printf("Here\n");
    value = x[0]*psisql(x[0],x[1])*sigma_gbs(x[0]);
    return value;
}

extern double uif_gbs_s (int *dim, double *x) {

    /*  r =  x[0]    
        z =  x[1]  */

    double value;
    value = x[0]*psisqs(x[0],x[1])*sigma_gbs(x[0]);
    return value;
}


extern double uif_gbs_c (int *dim, double *x) {

    /*  r =  x[0]    
        z =  x[1]  */

    double value;
    value = x[0]*psisqc(x[0],x[1])*sigma_gbw(x[0]);
    //value = x[0]*psisqc(x[0],x[1])*sigma_gbs(x[0]);
    return value;
}

extern double uif_gbs_b (int *dim, double *x) {

    /*  r =  x[0]    
        z =  x[1]  */

    double value;
    value = x[0]*psisqb(x[0],x[1])*sigma_gbs(x[0]);
    return value;
}


/*******************************************************************************
* The integrand of DIS GBW cross section
*******************************************************************************/
extern double uif_gbw (int *dim, double *x) {

    /*  r =  x[0]    
        z =  x[1]  
     */

    double value;

    switch (flavour) {
    case 0: /* Light  only */
        value = x[0]*psisql(x[0],x[1])*sigma_gbw(x[0]);
        break;
    case 1: /* Light + charm */
        value = x[0]*psisqlc(x[0],x[1])*sigma_gbw(x[0]);
        break;
    case 2:
        break;
    default:
        printf("Error: wrong 'flavour' value");
        break;
    }

   return value;
}

extern double uif_gbw_l (int *dim, double *x) {

    /*  r =  x[0]    
        z =  x[1]  */

    double value;
    value = x[0]*psisql(x[0],x[1])*sigma_gbw(x[0]);
    return value;
}

extern double uif_gbw_s (int *dim, double *x) {

    /*  r =  x[0]    
        z =  x[1]  */

    double value;
    value = x[0]*psisqs(x[0],x[1])*sigma_gbw(x[0]);
    return value;
}


extern double uif_gbw_c (int *dim, double *x) {

    /*  r =  x[0]    
        z =  x[1]  */

    double value;
    value = x[0]*psisqc(x[0],x[1])*sigma_gbw(x[0]);
    return value;
}

extern double uif_gbw_b (int *dim, double *x) {

    /*  r =  x[0]    
        z =  x[1]  */

    double value;
    value = x[0]*psisqb(x[0],x[1])*sigma_gbw(x[0]);
    return value;
}


/*******************************************************************************
* The integrand of DIS BGK cross section
*******************************************************************************/
extern double uif_bgk (int *dim, double *x) {

    /*  r =  x[0]    
        z =  x[1]  
     */

    double value,value1,value2;

    switch (flavour) {
    case 0: /* Light  only */
	switch (cheb_approx) {
	    case 0: /* Exact value - without Chebyshev approximation */
	    value = x[0]*psisql(x[0],x[1])*sigma_bgk(x[0]);
	    break;
	    case 1: /* With Chebyshev approximation */
	    value = x[0]*psisql(x[0],x[1])*sigma_bgk_cheb(x[0]);
	    break;
        } 
	break;
    case 1: /* Light + charm */
        //value = x[0]*psisqlc(x[0],x[1])*sigma_bgk_cheb(x[0]);
        light_charm  = 0;
        value1 = x[0]*psisql(x[0],x[1])*sigma_bgk_cheb(x[0]);
        light_charm  = 1;
        value2 = x[0]*psisqc(x[0],x[1])*sigma_bgk_cheb(x[0]);
	value=value1+value2;
        break;
    case 2: /* Charm only */
        value = x[0]*psisqc(x[0],x[1])*sigma_bgk_cheb(x[0]);
        break;
    default:
        printf("Error: wrong 'flavour' value");
        break;
    }

  // return value1 + value2;
   return value;
}

/*******************************************************************************
* The integrand of DIS BGK cross section  light quarks contribution
*******************************************************************************/
extern double uif_bgk_l (int *dim, double *x) {

    /*  r =  x[0]    
        z =  x[1]  
     */
    double value;
    value = x[0]*psisql(x[0],x[1])*sigma_bgk_cheb(x[0]);
    return value;
}

/*******************************************************************************
* The integrand of DIS BGK cross section  light quarks contribution
*******************************************************************************/
extern double uif_bgk_s (int *dim, double *x) {

    /*  r =  x[0]    
        z =  x[1]  
     */
    double value;
    value = x[0]*psisqs(x[0],x[1])*sigma_bgk_cheb(x[0]);
    return value;
}

/*******************************************************************************
* The integrand of DIS BGK cross section
*******************************************************************************/
extern double uif_bgk_c (int *dim, double *x) {

    /*  r =  x[0]    
        z =  x[1]  
     */
    double value;
    value = x[0]*psisqc(x[0],x[1])*sigma_bgk_cheb(x[0]);
    return value;
}

/*******************************************************************************
* The integrand of DIS BGK cross section
*******************************************************************************/
extern double uif_bgk_b (int *dim, double *x) {

    /*  r =  x[0]    
        z =  x[1]  
     */
    double value;
    value = x[0]*psisqb(x[0],x[1])*sigma_bgk_cheb(x[0]);
    return value;
}

/*******************************************************************************
* Unintegrated gluon distribution alpha*fpdf(x,l^2) integrand                  *
*******************************************************************************/
double alpha_fpdf_int (double *r) {

    double x_l, l2_l;
    double arg_j0;
    double norm_const, norm_units;
    double value;
    
    x_l  = arg_unint_x;
    l2_l = arg_unint_l2;

    arg_j0 = sqrt(l2_l)*r[0];
    
    xmod = x_l;

    norm_units = 1.0/0.389379292;            /* Units mb -> GeV{-2}*/
    norm_const = 3.0/(8.0*PI*PI)*l2_l*l2_l;

    switch (model) {
    case 0:
    value = norm_const*norm_units
            *r[0]*dbesj0_(&arg_j0)*(sigma_0-sigma_gbw(r[0]));
    break;
    case 1:
    value = norm_const*norm_units
            *r[0]*dbesj0_(&arg_j0)*(sigma_0-sigma_bgk_cheb(r[0]));
    break;
    case 2:
    value = norm_const*norm_units
            *r[0]*dbesj0_(&arg_j0)*(sigma_0-sigma_gbw(r[0]));
    break;
    }

    return value;
}

/*******************************************************************************
* Unintegrated gluon distribution alpha*fpdf(x,l^2)                            *
*******************************************************************************/
double alpha_fpdf (double x, double l2) {

    double rmi = 1.0e-05;
    double rma = 1.0e+04;
    double eps = 1.0e-09;

    extern double simps_();
    double dum1,dum2,dum3;
    double int_result;

    /* Copy x and l^2 to global variable */
    arg_unint_x = x;       
    arg_unint_l2 = l2;

   simps_(&rmi,&rma,&eps,&eps,&eps,
          &alpha_fpdf_int,&dum1,&int_result,&dum2,&dum3);
   //printf("%f  %f  %.10f\n",rmi, rma, int_result);

   return int_result;
}


/*******************************************************************************
* DIS cross section
*******************************************************************************/
double sigma_x (double X, double Q, double Y, double *par) {

    /* Integration limits */
    double A[2] = {amin,0};/* Lower limits of integration  {r,z} */ 
    double B[2] = {Rmax,0.5};       /* Upper limits of integration  {r,z} */ 


    /* dadmul integration variables */
    int    dim = 2;            /* Dimension of integral */
    int    minpts = 1.7e+01;
    int    maxpts = 1.0e+05;
    //int    maxpts = 1.0e+04;
    double eps = sigmaEPS;     /* Relative accuracy */
    //int    iwk = 2100;
    int    iwk = 21000;
    double wk[iwk];
    double result;             /* The result of integration */
    double relerr;             /* The relative accuracy of result */ 
    int    nfnevl;             /* The number of function eval. performed */
    int    ifail;              /* If '0' normal exit relerr < eps */
    //double temp;

    /* Set the global variabels according to sigma() arguments */    
    //xmod = X;
    xmod = X*(1+4*m_fsq/Q);
    xmod1 = X*(1+4*m_ch/Q);
    Q2   = Q;
     y   = Y;

    switch (model) {
        case 0:
        /* Set model parameters according to sigma() argument 'par' */    
        sigma_0 = par[0];
        lambda  = par[1];
        x_0     = par[2];

        /* Perform integration */
        //dadmul_(&uif_gbw, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
        //    &result, &relerr, &nfnevl, &ifail);

        B[0] = Rmax_simps2d;
        result = simps2d(A,B,eps_simps2d, uif_gbw);
        break;
       
        case 1:
        /* Set model parameters according to sigma() argument 'par' */    
        sigma_0  = par[0]; 
        A_g      = par[1];
        lambda_g = par[2];
        C        = par[3]; 
        mu02     = par[4];

       
        /* Perform integration */

        switch  (uif_int) {
           /* dadmul method */
           case 0:
           /*
           B[0] = 5.0;
           A[0] = 0.0;
           A[1] = 0.0;
           dadmul_(&uif_bgk, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
                   &result, &relerr, &nfnevl, &ifail);
           temp = result;
           */
           B[0] = Rmax_dadmul;
           //A[0] = 5.0;
           A[1] = 0.0;
           dadmul_(&uif_bgk, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
                   &result, &relerr, &nfnevl, &ifail);

          //result = result + temp;
               /* Print error message when ifail>0 */
               if(ifail>0)
                    printf("Error: ifail = %d \n", ifail); 
           break;

           case 1:
           /* simps2d method */
           B[0] = Rmax_simps2d;
           result = simps2d(A,B,eps_simps2d, uif_bgk);
    //printf("%e  %e  %e   %e   %e\n",sigma_0,A_g,lambda_g,C,mu02);
           break;
        }
        break;

        case 2:
        /* Set model parameters according to sigma() argument 'par' */    
        sigma_0 = par[0];
        lambda  = par[1];
        x_0     = par[2];
        C       = par[3];
	//mu02    = par[4];

	r_max	=par[4];
	//bmax	=par[4];
	g1      = par[5];
	g2      = par[6];

        /* Perform integration */
        //dadmul_(&uif_gbw, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
        //    &result, &relerr, &nfnevl, &ifail);

        B[0] = Rmax_simps2d;
        result = simps2d(A,B,eps_simps2d, uif_gbs);
        break;
    }         

    
    //printf("%2.10e  %2.8e  %d   %d \n",result,relerr,ifail,nfnevl);

    return result;
}

/*******************************************************************************
* DIS cross section
*******************************************************************************/
double sigma_l3 (double X, double Q, double Y, double *par) {

    /* Integration limits */
    double A[2] = {amin,0.0};/* Lower limits of integration  {r,z} */ 
    double B[2] = {Rmax,0.5};       /* Upper limits of integration  {r,z} */ 


    /* dadmul integration variables */
    int    dim = 2;            /* Dimension of integral */
    int    minpts = 1.7e+01;
    int    maxpts = 1.0e+05;
    //int    maxpts = 1.0e+04;
    double eps = sigmaEPS;     /* Relative accuracy */
    //int    iwk = 2100;
    int    iwk = 21000;
    double wk[iwk];
    double result;             /* The result of integration */
    double relerr;             /* The relative accuracy of result */ 
    int    nfnevl;             /* The number of function eval. performed */
    int    ifail;              /* If '0' normal exit relerr < eps */
    //double temp;

    /* Set the global variabels according to sigma() arguments */    
    //xmod = X;
    switch (xbj_mod){
    	case 0:
	xmod = X;
    	break;  
    	case 1:
	xmod = X*(1+4*m_fsq/Q);
    	break;  
    }

    Q2  = Q;
    y   = Y;

    switch (model) {
        case 0:
        /* Set model parameters according to sigma() argument 'par' */    
        sigma_0 = par[0];
        lambda  = par[1];
        x_0     = par[2];

        /* Perform integration */
        //dadmul_(&uif_gbw, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
        //    &result, &relerr, &nfnevl, &ifail);

        B[0] = Rmax_simps2d;
        result = simps2d(A,B,eps_simps2d, uif_gbw_l);
        break;
       
        case 1:
        /* Set model parameters according to sigma() argument 'par' */    
        sigma_0  = par[0]; 
        A_g      = par[1];
        lambda_g = par[2];
        C        = par[3]; 
        mu02     = par[4];

       
        /* Perform integration */

        switch  (uif_int) {
           /* dadmul method */
           case 0:
           /*
           B[0] = 5.0;
           A[0] = 0.0;
           A[1] = 0.0;
           dadmul_(&uif_bgk, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
                   &result, &relerr, &nfnevl, &ifail);
           temp = result;
           */
           B[0] = Rmax_dadmul;
           //A[0] = 5.0;
           A[1] = 0.0;
           dadmul_(&uif_bgk, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
                   &result, &relerr, &nfnevl, &ifail);

          //result = result + temp;
               /* Print error message when ifail>0 */
               if(ifail>0)
                    printf("Error: ifail = %d \n", ifail); 
           break;

           case 1:
           /* simps2d method */
           B[0] = Rmax_simps2d;
           result = simps2d(A,B,eps_simps2d, uif_bgk_l);
    //printf("%e  %e  %e   %e   %e\n",sigma_0,A_g,lambda_g,C,mu02);
           break;
        }
        break;
        case 2:
        /* Set model parameters according to sigma() argument 'par' */    
        sigma_0 = par[0];
        lambda  = par[1];
        x_0     = par[2];
        C       = par[3]; 
        //mu02    = par[4];
<<<<<<< HEAD
	r_max	= par[4];
=======
	bmax    = par[4];
>>>>>>> 5c56e8f61c56da62d800909f013271b2746be3fc
	g1      = par[5];
	g2      = par[6];

        /* Perform integration */
        //dadmul_(&uif_gbw, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
        //    &result, &relerr, &nfnevl, &ifail);

        B[0] = Rmax_simps2d;
        result = simps2d(A,B,eps_simps2d, uif_gbs_l);
        break;
    }         

    
    //printf("%2.10e  %2.8e  %d   %d \n",result,relerr,ifail,nfnevl);

    //printf("%f\n", result);
    return 6.0/5.0*result;
    //return 0.0;
}



/*******************************************************************************
* DIS cross section
*******************************************************************************/
double sigma_l (double X, double Q, double Y, double *par) {

    //printf("sigma_l");
    /* Integration limits */
    double A[2] = {amin,0};/* Lower limits of integration  {r,z} */ 
    double B[2] = {Rmax,0.5};       /* Upper limits of integration  {r,z} */ 


    /* dadmul integration variables */
    int    dim = 2;            /* Dimension of integral */
    int    minpts = 1.7e+01;
    int    maxpts = 1.0e+05;
    //int    maxpts = 1.0e+04;
    double eps = sigmaEPS;     /* Relative accuracy */
    //int    iwk = 2100;
    int    iwk = 21000;
    double wk[iwk];
    double result;             /* The result of integration */
    double relerr;             /* The relative accuracy of result */ 
    int    nfnevl;             /* The number of function eval. performed */
    int    ifail;              /* If '0' normal exit relerr < eps */
    //double temp;

    /* Set the global variabels according to sigma() arguments */    
    //xmod = X;
    switch (xbj_mod){
    	case 0:
	xmod = X;
    	break;  
    	case 1:
	xmod = X*(1+4*m_fsq/Q);
    	break;  
    }

    Q2  = Q;
    y   = Y;

    switch (model) {
        case 0:
        /* Set model parameters according to sigma() argument 'par' */    
        sigma_0 = par[0];
        lambda  = par[1];
        x_0     = par[2];

        /* Perform integration */
        //dadmul_(&uif_gbw, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
        //    &result, &relerr, &nfnevl, &ifail);

        B[0] = Rmax_simps2d;
        result = simps2d(A,B,eps_simps2d, uif_gbw_l);
        break;
       
        case 1:
        /* Set model parameters according to sigma() argument 'par' */    
        sigma_0  = par[0]; 
        A_g      = par[1];
        lambda_g = par[2];
        C        = par[3]; 
        mu02     = par[4];

       
        /* Perform integration */

        switch  (uif_int) {
           /* dadmul method */
           case 0:
           /*
           B[0] = 5.0;
           A[0] = 0.0;
           A[1] = 0.0;
           dadmul_(&uif_bgk, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
                   &result, &relerr, &nfnevl, &ifail);
           temp = result;
           */
           B[0] = Rmax_dadmul;
           //A[0] = 5.0;
           A[1] = 0.0;
           dadmul_(&uif_bgk, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
                   &result, &relerr, &nfnevl, &ifail);

          //result = result + temp;
               /* Print error message when ifail>0 */
               if(ifail>0)
                    printf("Error: ifail = %d \n", ifail); 
           break;

           case 1:
           /* simps2d method */
           B[0] = Rmax_simps2d;
           result = simps2d(A,B,eps_simps2d, uif_bgk_l);
    //printf("%e  %e  %e   %e   %e\n",sigma_0,A_g,lambda_g,C,mu02);
           break;
        }
        break;
        case 2:
        /* Set model parameters according to sigma() argument 'par' */    
        sigma_0 = par[0];
        lambda  = par[1];
        x_0     = par[2];
        C       = par[3]; 
        //mu02    = par[4];
	r_max	= par[4];
	//bmax    = par[4];

	g1      = par[5];
       	g2      = par[6];	

        /* Perform integration */
        //dadmul_(&uif_gbw, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
        //    &result, &relerr, &nfnevl, &ifail);

        B[0] = Rmax_simps2d;
        result = simps2d(A,B,eps_simps2d, uif_gbs_l);
        break;
       
    }         

    
    //printf("%2.10e  %2.8e  %d   %d \n",result,relerr,ifail,nfnevl);

    return result;
    //return 0.0;
}


/*******************************************************************************
* DIS cross section
*******************************************************************************/
double sigma_s (double X, double Q, double Y, double *par) {

    /* Integration limits */
    double A[2] = {amin,0};/* Lower limits of integration  {r,z} */ 
    double B[2] = {Rmax,0.5};       /* Upper limits of integration  {r,z} */ 


    /* dadmul integration variables */
    int    dim = 2;            /* Dimension of integral */
    int    minpts = 1.7e+01;
    int    maxpts = 1.0e+05;
    //int    maxpts = 1.0e+04;
    double eps = sigmaEPS;     /* Relative accuracy */
    //int    iwk = 2100;
    int    iwk = 21000;
    double wk[iwk];
    double result;             /* The result of integration */
    double relerr;             /* The relative accuracy of result */ 
    int    nfnevl;             /* The number of function eval. performed */
    int    ifail;              /* If '0' normal exit relerr < eps */
    //double temp;

    /* Set the global variabels according to sigma() arguments */    
    //xmod = X;
    switch (xbj_mod){
    	case 0:
	xmod = X;
    	break;  
    	case 1:
	xmod = X*(1+4*m_s/Q);
    	break;  
    }

    Q2  = Q;
    y   = Y;

    switch (model) {
        case 0:
        /* Set model parameters according to sigma() argument 'par' */    
        sigma_0 = par[0];
        lambda  = par[1];
        x_0     = par[2];

        /* Perform integration */
        //dadmul_(&uif_gbw, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
        //    &result, &relerr, &nfnevl, &ifail);

        B[0] = Rmax_simps2d;
        result = simps2d(A,B,eps_simps2d, uif_gbw_s);
        break;
       
        case 1:
        /* Set model parameters according to sigma() argument 'par' */    
        sigma_0  = par[0]; 
        A_g      = par[1];
        lambda_g = par[2];
        C        = par[3]; 
        mu02     = par[4];

       
        /* Perform integration */

        switch  (uif_int) {
           /* dadmul method */
           case 0:
           /*
           B[0] = 5.0;
           A[0] = 0.0;
           A[1] = 0.0;
           dadmul_(&uif_bgk, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
                   &result, &relerr, &nfnevl, &ifail);
           temp = result;
           */
           B[0] = Rmax_dadmul;
           //A[0] = 5.0;
           A[1] = 0.0;
           dadmul_(&uif_bgk, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
                   &result, &relerr, &nfnevl, &ifail);

          //result = result + temp;
               /* Print error message when ifail>0 */
               if(ifail>0)
                    printf("Error: ifail = %d \n", ifail); 
           break;

           case 1:
           /* simps2d method */
           B[0] = Rmax_simps2d;
           result = simps2d(A,B,eps_simps2d, uif_bgk_s);
    //printf("%e  %e  %e   %e   %e\n",sigma_0,A_g,lambda_g,C,mu02);
           break;
        }
        break;
        case 2:
        /* Set model parameters according to sigma() argument 'par' */    
        sigma_0 = par[0];
        lambda  = par[1];
        x_0     = par[2];
        C       = par[3]; 
        //mu02    = par[4];
        //mu02	= par[4];
	r_max    = par[4];
        g1      = par[5]; 
	g2	= par[6];
 	/* Perform integration */
        //dadmul_(&uif_gbw, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
        //    &result, &relerr, &nfnevl, &ifail);

        B[0] = Rmax_simps2d;
        result = simps2d(A,B,eps_simps2d, uif_gbs_s);
        break;
       
    }         

    
    //printf("%2.10e  %2.8e  %d   %d \n",result,relerr,ifail,nfnevl);

    return result;
}



/*******************************************************************************
* DIS cross section
*******************************************************************************/
double sigma_c (double X, double Q, double Y, double *par) {

    /* Integration limits */
    double A[2] = {amin,0.0};/* Lower limits of integration  {r,z} */ 
    double B[2] = {Rmax,0.5};       /* Upper limits of integration  {r,z} */ 


    /* dadmul integration variables */
    int    dim = 2;            /* Dimension of integral */
    int    minpts = 1.7e+01;
    int    maxpts = 1.0e+05;
    //int    maxpts = 1.0e+04;
    double eps = sigmaEPS;     /* Relative accuracy */
    //int    iwk = 2100;
    int    iwk = 21000;
    double wk[iwk];
    double result;             /* The result of integration */
    double relerr;             /* The relative accuracy of result */ 
    int    nfnevl;             /* The number of function eval. performed */
    int    ifail;              /* If '0' normal exit relerr < eps */
    //double temp;

    /* Set the global variabels according to sigma() arguments */    
    //xmod = X;
    switch (xbj_mod){
    	case 0:
	xmod = X;
    	break;  
    	case 1:
        if (photo ==1)    
	    xmod = 4*m_ch/(209.0*209.0);
        else
	    xmod = X*(1+4*m_ch/Q);
    	break;  
    }

    Q2   = Q;
     y   = Y;

    //printf("%f\n",Q);

    switch (model) {
        case 0:
        /* Set model parameters according to sigma() argument 'par' */    
        sigma_0 = par[0];
        lambda  = par[1];
        x_0     = par[2];

        /* Perform integration */
        //dadmul_(&uif_gbw, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
        //    &result, &relerr, &nfnevl, &ifail);

        B[0] = Rmax_simps2d;
        result = simps2d(A,B,eps_simps2d, uif_gbw_c);
        break;
       
        case 1:
        /* Set model parameters according to sigma() argument 'par' */    
        sigma_0  = par[0]; 
        A_g      = par[1];
        lambda_g = par[2];
        C        = par[3]; 
        mu02     = par[4];

       
        /* Perform integration */

        switch  (uif_int) {
           /* dadmul method */
           case 0:
           /*
           B[0] = 5.0;
           A[0] = 0.0;
           A[1] = 0.0;
           dadmul_(&uif_bgk, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
                   &result, &relerr, &nfnevl, &ifail);
           temp = result;
           */
           B[0] = Rmax_dadmul;
           //A[0] = 5.0;
           A[1] = 0.0;
           dadmul_(&uif_bgk, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
                   &result, &relerr, &nfnevl, &ifail);

          //result = result + temp;
               /* Print error message when ifail>0 */
               if(ifail>0)
                    printf("Error: ifail = %d \n", ifail); 
           break;

           case 1:
           /* simps2d method */
           B[0] = Rmax_simps2d;
           result = simps2d(A,B,eps_simps2d, uif_bgk_c);
    //printf("%e  %e  %e   %e   %e\n",sigma_0,A_g,lambda_g,C,mu02);
           break;
        }
        break;

        case 2:
        /* Set model parameters according to sigma() argument 'par' */    
        sigma_0 = par[0];
        lambda  = par[1];
        x_0     = par[2];
        C       = par[3]; 
        //mu02    = par[4];
	r_max	= par[4];
	//bmax	= par[4];
	g1      = par[5];
	g2      = par[6];

        /* Perform integration */
        //dadmul_(&uif_gbw, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
        //    &result, &relerr, &nfnevl, &ifail);

        B[0] = Rmax_simps2d;
        result = simps2d(A,B,eps_simps2d, uif_gbs_c);
        break;
    }         

    
    //printf("%2.10e  %2.8e  %d   %d \n",result,relerr,ifail,nfnevl);

    return result;
}

/*******************************************************************************
* DIS cross section
*******************************************************************************/
double sigma_b (double X, double Q, double Y, double *par) {

    /* Integration limits */
    double A[2] = {amin,0.0};/* Lower limits of integration  {r,z} */ 
    double B[2] = {Rmax,0.5};       /* Upper limits of integration  {r,z} */ 


    /* dadmul integration variables */
    int    dim = 2;            /* Dimension of integral */
    int    minpts = 1.7e+01;
    int    maxpts = 1.0e+05;
    //int    maxpts = 1.0e+04;
    double eps = sigmaEPS;     /* Relative accuracy */
    //int    iwk = 2100;
    int    iwk = 21000;
    double wk[iwk];
    double result;             /* The result of integration */
    double relerr;             /* The relative accuracy of result */ 
    int    nfnevl;             /* The number of function eval. performed */
    int    ifail;              /* If '0' normal exit relerr < eps */
    //double temp;

    /* Set the global variabels according to sigma() arguments */    
    //xmod = X;
    switch (xbj_mod){
    	case 0:
	xmod = X;
    	break;  
    	case 1:
	if (photo ==1)    
	    xmod = 4*m_b/(209*209);
        else
	    xmod = X*(1+4*m_b/Q);
    	break;  
    }

    Q2   = Q;
     y   = Y;

    //printf("%f\n",Q);

    switch (model) {
        case 0:
        /* Set model parameters according to sigma() argument 'par' */    
        sigma_0 = par[0];
        lambda  = par[1];
        x_0     = par[2];

        /* Perform integration */
        //dadmul_(&uif_gbw, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
        //    &result, &relerr, &nfnevl, &ifail);

        B[0] = Rmax_simps2d;
        result = simps2d(A,B,eps_simps2d, uif_gbw_b);
        break;
       
        case 1:
        /* Set model parameters according to sigma() argument 'par' */    
        sigma_0  = par[0]; 
        A_g      = par[1];
        lambda_g = par[2];
        C        = par[3]; 
        mu02     = par[4];

       
        /* Perform integration */

        switch  (uif_int) {
           /* dadmul method */
           case 0:
           /*
           B[0] = 5.0;
           A[0] = 0.0;
           A[1] = 0.0;
           dadmul_(&uif_bgk, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
                   &result, &relerr, &nfnevl, &ifail);
           temp = result;
           */
           B[0] = Rmax_dadmul;
           //A[0] = 5.0;
           A[1] = 0.0;
           dadmul_(&uif_bgk, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
                   &result, &relerr, &nfnevl, &ifail);

          //result = result + temp;
               /* Print error message when ifail>0 */
               if(ifail>0)
                    printf("Error: ifail = %d \n", ifail); 
           break;

           case 1:
           /* simps2d method */
           B[0] = Rmax_simps2d;
           result = simps2d(A,B,eps_simps2d, uif_bgk_b);
    //printf("%e  %e  %e   %e   %e\n",sigma_0,A_g,lambda_g,C,mu02);
           break;
        }
        break;

        case 2:
        /* Set model parameters according to sigma() argument 'par' */    
        sigma_0 = par[0];
        lambda  = par[1];
        x_0     = par[2];
        C       = par[3]; 
        //mu02    = par[4];
	r_max 	=par[4];
	//bmax    = par[4];
	g1      = par[5]; 
	g2	= par[6];

        /* Perform integration */
        //dadmul_(&uif_gbw, &dim, &A, &B, &minpts, &maxpts, &eps, &wk, &iwk,
        //    &result, &relerr, &nfnevl, &ifail);

        B[0] = Rmax_simps2d;
        result = simps2d(A,B,eps_simps2d, uif_gbs_b);
        break;
    }         

    
    //printf("%2.10e  %2.8e  %d   %d \n",result,relerr,ifail,nfnevl);

    return result;
}


/*******************************************************************************
* The function which returns maximal of two arguments
*******************************************************************************/
double sigma (double X, double Q, double Y, double *par) {

    //printf("%f %f\n", Q2, Q);
    double value;
    switch (flavour) {
        case 0:
        //value = sigma_l (X,Q,Y,par);
        value = sigma_l (X,Q,Y,par) + sigma_s (X,Q,Y,par);
        break;
        case 1:
        if (fl_beauty==0) 
	    value = sigma_l3 (X,Q,Y,par) + sigma_c (X,Q,Y,par);
        else if (fl_beauty==1) {
	    value = sigma_l3 (X,Q,Y,par) + sigma_c (X,Q,Y,par)
	           + sigma_b (X,Q,Y,par); }
        break;
        case 2:
        value = sigma_c (X,Q,Y,par); 
        break;
        case 3:
        value = sigma_b (X,Q,Y,par); 
        break;
    }
    //printf("%f\n", value);
    return value;
}

/*******************************************************************************
* The function which returns maximal of two arguments
*******************************************************************************/
double maximum (double arg1, double arg2) {

    double value;
    if(arg1>arg2)
	value = arg1;
    else if(arg1<=arg2)
	value = arg2;
    return value;
}

/*******************************************************************************
* The function counting number of data
*******************************************************************************/
int count_data (char * file_name) {

    FILE *datafile; 
    double dum;
    double stat,sysup,sysdown;

     int good_data = 0;
     int i;
 
    //double q_down = 0.0;
    //double q_up   = 4000;
    //double x_up   = 0.01;
   

    //printf("%s\n",file_name);
    datafile = fopen(file_name,"r");    
    /* Read experimental data from file */
    if(strcmp(file_name,"data/h1zeus_05.dat")==0) {
	for(i=0;i<331;i++) {
	    fscanf(datafile,
		   " %lE %lE %lE %lE %lE", 
		   &xbdata[i], &ydata[i], &q2data[i], &csdata[i], &erdata[i]);

	    if ((xbdata[i]<=x_up)&&(q2data[i]>=q_down)&&(q2data[i]<=q_up)) {
		good_data++;
	    }
	}
    }
    else if(strcmp(file_name,"h1low2001.dat")==0) {
	for(i=0;i<62;i++) {
	    fscanf(datafile,
		   "   %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE", 
		   &q2data[i],&xbdata[i],&ydata[i],&dum,&dum, 
                   &csdata[i], &erdata[i],&dum,&dum,&dum);
            erdata[i]= erdata[i]*csdata[i]/100;

	    if ((xbdata[i]<=x_up)&&(q2data[i]>=q_down)&&(q2data[i]<=q_up)) {
		good_data++;
	    }
        }

    } 
    else if(strcmp(file_name,"h12001.dat")==0) {
	for(i=0;i<71;i++) {
	    fscanf(datafile,
		   "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
		   &q2data[i],&xbdata[i],&ydata[i],&dum,&dum, 
                   &csdata[i], &erdata[i],&dum,&dum,&dum);
            erdata[i]= erdata[i]*csdata[i]/100;

	    if ((xbdata[i]<=x_up)&&(q2data[i]>=q_down)&&(q2data[i]<=q_up)) {
		good_data++;
	    }
            //printf("%f %f %f %f\n", 
	    //	   q2data[i],xbdata[i],ydata[i],csdata[i]);

        }

    }
    else if(strcmp(file_name,"bpt97.dat")==0) {
	for(i=0;i<70;i++) {
	    fscanf(datafile,
		   "%lf %lE %lf %lf %lf %lf %lf", 
		   &q2data[i],&xbdata[i],&ydata[i],&csdata[i],
                   &stat,&sysup,&sysdown); 
            dum = maximum(sysup,sysdown);
            erdata[i] = sqrt(stat*stat+dum*dum);

	    if ((xbdata[i]<=x_up)&&(q2data[i]>=q_down)&&(q2data[i]<=q_up)) {
		good_data++;
	    }
        }

    }
    else if(strcmp(file_name,"zeus2001.dat")==0) {
        //int k = 1;
	for(i=0;i<242;i++) {
	    fscanf(datafile,
		   "    %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
		   &q2data[i],&xbdata[i],&csdata[i],&dum, 
                   &sysup,&sysdown,&stat,
                   &dum,&dum,&dum,&dum,&dum,&dum,&dum,&dum,&dum,&dum,&dum); 
            stat  = stat*csdata[i]/100;
            sysup = sysup*csdata[i]/100;
            erdata[i] = sqrt(stat*stat+sysup*sysup);

	    if ((xbdata[i]<=x_up)&&(q2data[i]>=q_down)&&(q2data[i]<=q_up)) {
		good_data++;
/*
            printf("%3d  %e %e %e %e\n", 
		   k++, q2data[i],xbdata[i],csdata[i],erdata[i]);
*/
                   //stat,sysup,sysdown); 
	    }
        }
    }
    else if(strcmp(file_name,"data/hera_tot.dat")==0) {

	for(i=0;i<597;i++) {
	    fscanf(datafile,
		   "%lE %lE %lE %lE %lE", 
		   &q2data[i],&xbdata[i],&wdata[i],&csdata[i], &erdata[i]); 
	    if ((xbdata[i]<=x_up)&&(q2data[i]>=q_down)&&(q2data[i]<=q_up)) {
		good_data++;
	    }
        }

    }
    fclose(datafile);

    //printf("%d\n", good_data);
    return good_data;

}


/*******************************************************************************
* The function reading data file form unofficial 'h1zeus_05.dat' 
*******************************************************************************/
void readdata_h1zeus05(void) {

    FILE *datafile; 
     int i;

    datafile = fopen("data/h1zeus_05.dat","r");    

    /* Read experimental data from file */
    for(i=0;i<NDATA;i++) {

        fscanf(datafile,
               " %lE %lE %lE %lE %lE", 
               &xbdata[i], &ydata[i], &q2data[i], &csdata[i], &erdata[i]);
        switch (flavour) {
            case 0:
	    xmoddata[i] =  xbdata[i]*(1+4*m_fsq/q2data[i]);
            break;
            case 1:
	    xmoddata[i] =  xbdata[i];
	    //xmoddata[i] =  xbdata[i]*(1+4*m_ch/q2data[i]);
            break;
        }

        /* Ignore data with x>0.01 or q2>400  by overwriting arrays' elemts */
        if ((xbdata[i]>1.0e-2)||(q2data[i]>4000))
            i--;
       
//    printf("%e %e %e %e %e \n", 
//           xbdata[i], ydata[i], q2data[i], csdata[i], erdata[i]);
    }
    fclose(datafile);
}

/*******************************************************************************
* The function reading the most recent official data 
*******************************************************************************/
void readdata_new (void) {

    FILE *datafile; 
    int i;

    double dum;
    double stat,sysup,sysdown;

    int n_new; 

    n_new    = count_data("data/hera_tot.dat");

    printf("n_new:   %d\n",n_new);
    datafile = fopen("data/hera_tot.dat","r");    

   double alpha =7.297e-3;
   double xmp0 = 0.93827;

    for(i=0;i<n_new;i++) {
	//printf("%d %d \n",i, n_new);
	    fscanf(datafile,
		   "%lE %lE %lE %lE %lE", 
		   &q2data[i],&xbdata[i],&wdata[i],&csdata[i],&erdata[i]); 
//	    printf(">>> %lE %lE %lE %lE %lE\n ", 
//		   q2data[i],xbdata[i],wdata[i],csdata[i],erdata[i]); 
//
            double fac = pow(q2data[i],2)*(1-xbdata[i])/
                        (4*pow(pi,2)*alpha*(q2data[i]+pow(2*xbdata[i]*xmp0,2)));
            double units = 10.0 * pow(197.3271,2);
            fac = fac /(units * 1.e-3);    //[mikrobarn^-1]

            csdata[i] = fac*csdata[i];
            erdata[i] = fac*erdata[i];

//   fac = fac /(units * 1.d-3)    ![mikrobarn^-1]

//   fac = q2**2*(1-xb)/(4*pi**2*alpha*(q2+(2*xb*xmp0)**2))  ![GeV^2]
//   units = 10.d0 * (197.3271d0)**2  ![GeV^2 * nb]     convert GeV^-2 to nb
//   fac = fac /(units * 1.d-3)    ![mikrobarn^-1]


	    switch (flavour) {
		case 0:
		xmoddata[i] =  xbdata[i];
		//xmoddata[i] =  xbdata[i]*(1+4*m_fsq/q2data[i]);
		break;

		case 1:
		xmoddata[i] =  xbdata[i];
		//xmoddata[i] =  xbdata[i]*(1+4*m_fsq/q2data[i]);
		//xmoddata[i] =  xbdata[i]*(1+4*m_ch/q2data[i]);
		break;
	    }

        if ((xbdata[i]>x_up)||(q2data[i]<q_down)||(q2data[i]>q_up))
            i--;
        //printf("%e  %e  %e\n",xbdata[i], xmoddata[i], csdata[i]);
    } 

    nf2data = n_new;

    fclose(datafile);

}

/*******************************************************************************
* The function reading the most recent official data 
*******************************************************************************/
void readdata (void) {

    FILE *datafile; 
    int i;

    double dum;
    double stat,sysup,sysdown;

    //int n_05; 
    int n_h1low; 
    int n_h1; 
    int n_bpt; 
    int n_zeus; 

    //n_05    = count_data("h1zeus_05.dat");
    n_h1low = count_data("h1low2001.dat");
    n_h1    = count_data("h12001.dat");
    n_bpt   = count_data("bpt97.dat");
    n_zeus  = count_data("zeus2001.dat");
/*
    printf("%d\n",n_05);
*/
    printf("h1_low: %d\n",n_h1low);
    printf("h1:     %d\n",n_h1);
    printf("n_bpt:  %d\n",n_bpt);
    printf("zeus:   %d\n",n_zeus);
    printf("sum: %d\n",n_h1low+n_h1+n_bpt+n_zeus);
    datafile = fopen("data/h1low2001.dat","r");    

    for(i=0;i<n_h1low;i++) {
	    fscanf(datafile,
		   "   %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE", 
		   &q2data[i], &xbdata[i], &ydata[i],&dum,&dum, 
                   &csdata[i], &erdata[i],&dum,&dum,&dum);
            erdata[i]= rescale*erdata[i]*csdata[i]/100;
            csdata[i]= rescale*csdata[i];
	    switch (flavour) {
		case 0:
		xmoddata[i] =  xbdata[i];
		//xmoddata[i] =  xbdata[i]*(1+4*m_fsq/q2data[i]);
		break;

		case 1:
		xmoddata[i] =  xbdata[i];
		//xmoddata[i] =  xbdata[i]*(1+4*m_fsq/q2data[i]);
		//xmoddata[i] =  xbdata[i]*(1+4*m_ch/q2data[i]);
		break;
	    }

        if ((xbdata[i]>x_up)||(q2data[i]<q_down)||(q2data[i]>q_up))
            i--;
        //printf("%e  %e  %e\n",xbdata[i], xmoddata[i], csdata[i]);
    } 

    fclose(datafile);

    datafile = fopen("data/h12001.dat","r");    
    //for(i=0;i<n_h1;i++) {
    for(i=n_h1low;i<(n_h1+n_h1low);i++) {
	    fscanf(datafile,
		   "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
		   &q2data[i],&xbdata[i],&ydata[i],&dum,&dum, 
                   &csdata[i], &erdata[i],&dum,&dum,&dum);
            erdata[i]= rescale*erdata[i]*csdata[i]/100;
            csdata[i]= rescale*csdata[i];
	    switch (flavour) {
		case 0:
		xmoddata[i] =  xbdata[i];
		//xmoddata[i] =  xbdata[i]*(1+4*m_fsq/q2data[i]);
		break;

		case 1:
		xmoddata[i] =  xbdata[i];
		//xmoddata[i] =  xbdata[i]*(1+4*m_fsq/q2data[i]);
		//xmoddata[i] =  xbdata[i]*(1+4*m_ch/q2data[i]);
		break;
	    }

        if ((xbdata[i]>x_up)||(q2data[i]<q_down)||(q2data[i]>q_up))
            i--;
    }

    fclose(datafile);

    datafile = fopen("data/bpt97.dat","r");    
    for(i=(n_h1+n_h1low);i<(n_h1+n_h1low+n_bpt);i++) {
	    fscanf(datafile,
		   "%lf %lE %lf %lf %lf %lf %lf", 
		   &q2data[i],&xbdata[i],&ydata[i],&csdata[i],
                   &stat,&sysup,&sysdown); 
            dum = maximum(sysup,sysdown);
            erdata[i] = sqrt(stat*stat+dum*dum);
	    switch (flavour) {
		case 0:
		xmoddata[i] =  xbdata[i];
		//xmoddata[i] =  xbdata[i]*(1+4*m_fsq/q2data[i]);
		break;

		case 1:
		xmoddata[i] =  xbdata[i];
		//xmoddata[i] =  xbdata[i]*(1+4*m_fsq/q2data[i]);
		//xmoddata[i] =  xbdata[i]*(1+4*m_ch/q2data[i]);
		break;
	    }

        if ((xbdata[i]>x_up)||(q2data[i]<q_down)||(q2data[i]>q_up))
            i--;

    }

    fclose(datafile);
    datafile = fopen("data/zeus2001.dat","r");    
    for(i=(n_h1+n_h1low+n_bpt);i<(n_h1+n_h1low+n_bpt+n_zeus);i++) {
	    fscanf(datafile,
		   "    %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
		   &q2data[i],&xbdata[i],&csdata[i],&dum, 
                   &sysup,&sysdown,&stat,
                   &dum,&dum,&dum,&dum,&dum,&dum,&dum,&dum,&dum,&dum,&dum); 
            stat  = stat*csdata[i]/100;
            sysup = sysup*csdata[i]/100;
            erdata[i] = sqrt(stat*stat+sysup*sysup);
	    switch (flavour) {
		case 0:
		xmoddata[i] =  xbdata[i];
		//xmoddata[i] =  xbdata[i]*(1+4*m_fsq/q2data[i]);
		break;

		case 1:
		xmoddata[i] =  xbdata[i];
		//xmoddata[i] =  xbdata[i]*(1+4*m_fsq/q2data[i]);
		//xmoddata[i] =  xbdata[i]*(1+4*m_ch/q2data[i]);
		break;
	    }



        if ((xbdata[i]>x_up)||(q2data[i]<q_down)||(q2data[i]>q_up))
            i--;
    }

    nf2data = n_h1+n_h1low+n_bpt+n_zeus;

    fclose(datafile);

/*
    printf("No.       q2         xb           xmod          cs           er\n"); 
    printf("\n") ;
    //for(i=0;i<5;i++) {
    for(i=0;i<(n_h1+n_h1low+n_bpt+n_zeus);i++) {
            printf("%3d %e %e %e %e %e\n", 
                   i+1 , q2data[i], xbdata[i], xmoddata[i],csdata[i],erdata[i]);
    }

*/
}

/*******************************************************************************
* The function calculating ydata[] instead of taking experimental values
*******************************************************************************/
void  calc_ydata (void) {

    int    i;
    double s;

    s = 4*27.6*920; /* in GeV^2 */

    for(i=0;i<NDATA;i++)
	ydata[i] = q2data[i]/(xmoddata[i]*s);
}


/*******************************************************************************
* The function plotting F_2 inclusive 
*******************************************************************************/
void graphdata (void) {

    FILE *fout; 
     int i,k;

    char foutname[5];

   // double pmts[5] = {sigma_0,A_g,lambda_g,C,mu02};

    for(i=0;i<25;i++) {

        sprintf(foutname,"%.1f",q2graph[i]);
	fout = fopen(strcat(foutname,".gdt"),"w");    
	for(k=0;k<NDATA;k++) {

	    if (q2data[k]==q2graph[i]){
		   fprintf(fout,"%e  %e  %e  %e\n",xbdata[k], xmoddata[k],
		   //printf("%e  %e \n",xmoddata[k],
                   csdata[k], sfcn0[k]);
                   //        sigma(xmoddata[k],q2data[k],ydata[k],pmts));
            }
	}
	fclose(fout);
    }
}

/*******************************************************************************
* The function reading charm F_2 data file
*******************************************************************************
void readdata_charm (void) {

    FILE *h1x_file; 
    FILE *h1Q_file; 
    FILE *zeusx_file; 
    FILE *zeus!_file; 
     int i;

    h1x_file = fopen("h1x.dat","r");    
    h1Q_file = fopen("h1Q.dat","r");    
    zeusx_file = fopen("zeusx.dat","r");    
    zeusQ_file = fopen("zeusQ.dat","r");    

    for(i=0;i<NDATA;i++) {

        fscanf(datafile,
               " %lE %lE %lE %lE %lE", 
               &xbdata[i], &ydata[i], &q2data[i], &csdata[i], &erdata[i]);
        switch (flavour) {
            case 0:
	    xmoddata[i] =  xbdata[i]*(1+4*m_fsq/q2data[i]);
            break;

            case 1:
	    xmoddata[i] =  xbdata[i];
	    //xmoddata[i] =  xbdata[i]*(1+4*m_ch/q2data[i]);
            break;
        }

        if ((xbdata[i]>1.0e-2)||(q2data[i]>4000))
            i--;
       
//    printf("%e %e %e %e %e \n", 
//           xbdata[i], ydata[i], q2data[i], csdata[i], erdata[i]);
    }
    

    fclose(datafile);
}
*/

/*******************************************************************************
* 
*******************************************************************************/
void update_graph (char *graph_file) {

FILE *fp;
FILE *fk;

double par_1;
double par_2;
double par_3;
double par_4;
double par_5;

int c;

fp=fopen(graph_file,"r");    
fk=fopen("graph.tmp","w");

fscanf(fp,"par_1 = %le\n", &par_1);
fscanf(fp,"par_2 = %le\n", &par_2);
fscanf(fp,"par_3 = %le\n", &par_3);
fscanf(fp,"par_4 = %le\n", &par_4);
fscanf(fp,"par_5 = %le\n", &par_5);
 
par_1 = sigma_0;
par_2 = A_g;
par_3 = lambda_g;
par_4 = C;
par_5 = mu02;

fprintf(fk,"par_1 = %le\n", par_1);
fprintf(fk,"par_2 = %le\n", par_2);
fprintf(fk,"par_3 = %le\n", par_3);
fprintf(fk,"par_4 = %le\n", par_4);
fprintf(fk,"par_5 = %le\n", par_5);
fprintf(fk,"\n");

while((c=getc(fp))!=EOF) {        

    fprintf(fk,"%c",c);
} 

fclose(fp);
fclose(fk);

fp=fopen(graph_file,"w");    
fk=fopen("graph.tmp","r");

while((c=getc(fk))!=EOF) {        

    fprintf(fp,"%c",c);
} 

fclose(fp);
fclose(fk);
}

/*******************************************************************************
* The function plotting dipol cross section 
*******************************************************************************/
void graph_f_2_charm (void) {

    FILE *fout; 
     int i,k;

    double q2charm[14] = {1.5,2.0,3.5,4.0,6.5,7.0,11.0,12.0,18.0,25.0,30,60,
                        130.0,500.0};
    double x_Bj, f2c,f2c_mod;
    double par_temp[5];

    double npoints; /* Number of points used to plot */
    ////////////////double r = 0.0;


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

    npoints = 12;
    flavour = 2;

    chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);

	fout = fopen("graph_f_2_charm.dat","w");    

	for(i=0;i<14;i++) {
            for(k=0;k<npoints;k++) {
                //xmod = pow(10,-(0.5*k+1)); /* x = 10^-1...10^-6 */
                x_Bj = pow(10,-(0.5*k+1)); /* x = 10^-1...10^-6 */
		//xmod = x_Bj*(1+4*m_ch/q2charm[i]);
	        //printf("%e %e %e\n", 
		f2c     = sigma(x_Bj, q2charm[i],0.0,par_temp);
		f2c_mod = pow((1-x_Bj),20)*f2c;
	        fprintf(fout,"%e %e %e  %e\n", 
                      q2charm[i], x_Bj, f2c,f2c_mod);
                      //q2charm[i], xmod, sigma(xmod, q2charm[i],0.0,par_temp));
                      //q2charm[i], x_Bj, sigma(xmod, q2charm[i],0.0,par_temp));
                      //q2charm[i], x_Bj, xmod);
            }
	    fprintf(fout,"\n"); 
        }

	fclose(fout);

    update_graph("f_2_charm.plt");
}

/*******************************************************************************
* The function plotting 
*******************************************************************************/
void graph_f_2_charm_9 (int model_graph) {

    FILE *fout; 
     int i,k;

    double q2charm[9] = {2.0,4.0,7.0,11.0,18.0,30,60,130.0,500.0};
    double x_Bj, f2c,f2c_mod;
    double par_temp[5];

    double npoints; /* Number of points used to plot */
    ////////////////double r = 0.0;

   
    model = model_graph;

    switch (model) {
    case 0:
	sigma_0 = 29.12; 
        lambda  = 0.277; 
        x_0     = 0.41e-04; 
/*
	sigma_0 = par_start[0]; 
        lambda  = par_start[1]; 
        x_0     = par_start[2]; 
*/

	par_temp[0] = sigma_0; 
        par_temp[1] = lambda; 
        par_temp[2] = x_0; 
        printf("GBW\n");
	fout = fopen("f2_charm_gbw.dat","w");    
    break;
    case 1:
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

        printf("BGK\n");
	fout = fopen("f2_charm_bgk.dat","w");    
	chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);
    break;
    }

    npoints = 12;
    flavour = 2;


	//fout = fopen("graph_f_2_charm.dat","w");    

	for(i=0;i<9;i++) {
            for(k=0;k<npoints;k++) {
                //xmod = pow(10,-(0.5*k+1)); /* x = 10^-1...10^-6 */
                x_Bj = pow(10,-(0.5*k+1)); /* x = 10^-1...10^-6 */
		//xmod = x_Bj*(1+4*m_ch/q2charm[i]);
	        //printf("%e %e %e\n", 
		f2c     = sigma(x_Bj, q2charm[i],0.0,par_temp);
		f2c_mod = pow((1-x_Bj),20)*f2c;
	        fprintf(fout,"%e %e %e  %e\n", 
                      q2charm[i], x_Bj, f2c,f2c_mod);
                      //q2charm[i], xmod, sigma(xmod, q2charm[i],0.0,par_temp));
                      //q2charm[i], x_Bj, sigma(xmod, q2charm[i],0.0,par_temp));
                      //q2charm[i], x_Bj, xmod);
            }
	    fprintf(fout,"\n"); 
        }

	fclose(fout);

    update_graph("f_2_charm.plt");
}

/*******************************************************************************
* The function plotting 
*******************************************************************************/
void graph_f2ch_comp (void) {

    FILE *fout; 
     int i,k;

    double q2charm[16] = 
    {1.5,2.0,2.5,3.5,5.0,6.5,8.5,12.0,15.0,20.0,25.0,35.0,45.0,60.0,90.0,120.0};
    double x_Bj, f2c[16], f2a[16];
    double par_temp[5];

    double npoints; /* Number of points used to plot */
   
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

    printf("BGK for comparision with DGLAP\n");

    fout = fopen("f2ch_comp.dat","w");    
    chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);

    npoints = 50;

    for(k=0;k<npoints;k++) {
	x_Bj = pow(10,-(4.0*k/(npoints-1)+1.0)); /* x = 10^-1...10^-6 */

	flavour = 1;
	for(i=0;i<16;i++) {
	    f2a[i]  = sigma(x_Bj, q2charm[i],0.0,par_temp);
            if(negglu==1) 
		f2a[i] = sigma(pow(10,-(4.0/(npoints-1)+1.0))*x_Bj, 
                               q2charm[i],0.0,par_temp);
	}
	flavour = 2;
	for(i=0;i<16;i++) {
	    f2c[i]  = sigma(x_Bj, q2charm[i],0.0,par_temp);
	}
	printf("%.3E   %E   %E   %E   %E\n",x_Bj,f2a[0],f2a[15],f2c[0],f2c[15]);
	fprintf(fout,
               "%E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E", x_Bj,
	       f2a[0],f2a[1],f2a[2],f2a[3],f2a[4],f2a[5],f2a[6],f2a[7],f2a[8],
	       f2a[9],f2a[10],f2a[11],f2a[12],f2a[13],f2a[14],f2a[15]);
	fprintf(fout,
               "%E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E\n",
	       f2c[0],f2c[1],f2c[2],f2c[3],f2c[4],f2c[5],f2c[6],f2c[7],f2c[8],
	       f2c[9],f2c[10],f2c[11],f2c[12],f2c[13],f2c[14],f2c[15]);
    }
    fclose(fout);
}
/*******************************************************************************
* The function plotting 
*******************************************************************************/
void graph_f2c9q2 (int model_graph) {

    FILE *fout; 
     int i,k;

    //double x_Bj[6] = {0.0002,0.0005,0.002,0.005,0.013,0.032};
    double x_Bj[6] = {0.032,0.013,0.005,0.002,0.0005,0.0002};
    double q2charm, f2c;
    double par_temp[5];

    double npoints; /* Number of points used to plot */
    ////////////////double r = 0.0;

   
    model = model_graph;

    switch (model) {
    case 0:
/*
	sigma_0 = 29.12; 
        lambda  = 0.277; 
        x_0     = 0.41e-04; 
*/

	sigma_0 = par_start[0]; 
        lambda  = par_start[1]; 
        x_0     = par_start[2]; 

	par_temp[0] = sigma_0; 
        par_temp[1] = lambda; 
        par_temp[2] = x_0; 
        printf("GBW\n");
	fout = fopen("f2_chgbwq2.dat","w");    
    break;
    case 1:
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

        printf("BGK\n");
	fout = fopen("f2_chbgkq2.dat","w");    
	chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);
    break;
    }

    npoints = 12;
    flavour = 2;


	for(i=0;i<6;i++) {
            for(k=0;k<npoints;k++) {
                q2charm = pow(10,(0.28*k)); /* */
		f2c = pow(4,i)*sigma(x_Bj[i], q2charm, 0.0,par_temp);
	        fprintf(fout,"%e  %e    %e\n",x_Bj[i], q2charm, f2c);
	        printf("%e  %e    %e\n",x_Bj[i], q2charm, f2c);
            }
	    fprintf(fout,"\n"); 
        }

	fclose(fout);


		//f2c     = sigma(0.0005, 25, 0.0,par_temp);
	        //printf("%e  %e    %e\n",x_Bj[i], q2charm, f2c);

    update_graph("f_2_charm.plt");
}

/*******************************************************************************
* The function plotting 
*******************************************************************************/
void graph_f2b9q2 (int model_graph) {

    FILE *fout; 
     int i,k;

    double x_Bj[6] = {0.032,0.013,0.005,0.002,0.0005,0.0002};
    double q2beauty, f2b;
    double par_temp[5];

    double npoints; /* Number of points used to plot */

   
    model = model_graph;

    switch (model) {
    case 0:
/*
	sigma_0 = 29.12; 
        lambda  = 0.277; 
        x_0     = 0.41e-04; 
*/

	sigma_0 = par_start[0]; 
        lambda  = par_start[1]; 
        x_0     = par_start[2]; 

	par_temp[0] = sigma_0; 
        par_temp[1] = lambda; 
        par_temp[2] = x_0; 
        printf("GBW\n");
	fout = fopen("f2_begbwq2.dat","w");    
    break;
    case 1:
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

        printf("BGK\n");
	fout = fopen("f2_bebgkq2.dat","w");    
	chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);
    break;
    }

    npoints = 12;
    flavour = 3;


	for(i=0;i<6;i++) {
            for(k=0;k<npoints;k++) {
                if((x_Bj[i]==0.032)&&model==1)
		    q2beauty = 9*pow(10,(0.30*k)); /* */
                else if((x_Bj[i]==0.013)&&model==1)
		    q2beauty = 5*pow(10,(0.30*k)); /* */
                else
		    q2beauty = pow(10,(0.28*k)); /* */
		f2b = pow(8,i)*sigma(x_Bj[i], q2beauty, 0.0,par_temp);
	        fprintf(fout,"%e  %e    %e\n",x_Bj[i], q2beauty, f2b);
	        printf("%e  %e    %e\n",x_Bj[i], q2beauty, f2b);
            }
	    fprintf(fout,"\n"); 
        }

	fclose(fout);

}


/*******************************************************************************
* The function plotting 
*******************************************************************************/
void graph_f2chf2 (void) {

    FILE *fout; 
     int i,k,m = 1;

    double q2charm[4] = {2.0,10.0,50.0,100.0};
    double x_Bj, f2c[4][3], f2a[4][3];
    double par_temp[5];

    double npoints; /* Number of points used to plot */
   

    model = 1;
    switch(model) {
    case 0:
	sigma_0 = par_start[0]; 
	lambda  = par_start[1]; 
	x_0     = par_start[2]; 

	par_temp[0] = sigma_0; 
	par_temp[1] = lambda; 
	par_temp[2] = x_0; 
        break;
    case 1:
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

        break;
    }

    fout = fopen("f2chf2.dat","w");    

    npoints = 10;

    for(k=0;k<npoints;k++) {
	x_Bj = pow(10,-(4.0*k/(npoints-1)+1.0)); /* x = 10^-1...10^-6 */

        /* 1st set of parameters */
	flavour = 1;
	par_temp[m] = bgk_parst[m]; 
	chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);
	//par_temp[m] = par_start[m]; 
	for(i=0;i<4;i++) {
	    f2a[i][0]  = sigma(x_Bj, q2charm[i],0.0,par_temp);
	}
	flavour = 2;
	for(i=0;i<4;i++) {
	    f2c[i][0]  = sigma(x_Bj, q2charm[i],0.0,par_temp);
	}

        /* 2nd set of parameters */
	flavour = 1;
	par_temp[m] = bgk_parst[m]/2.0; 
	chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);
	//par_temp[m] = par_start[m]/2.0; 
	for(i=0;i<4;i++) {
	    f2a[i][1]  = sigma(x_Bj, q2charm[i],0.0,par_temp);
	}
	flavour = 2;
	for(i=0;i<4;i++) {
	    f2c[i][1]  = sigma(x_Bj, q2charm[i],0.0,par_temp);
	}

	/* 3rd set of parameters */
	flavour = 1;
	par_temp[m] = bgk_parst[m]/10.0; 
	chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);
	//par_temp[m] = par_start[m]/10.0; 
	for(i=0;i<4;i++) {
	    f2a[i][2]  = sigma(x_Bj, q2charm[i],0.0,par_temp);
	}
	flavour = 2;
	for(i=0;i<4;i++) {
	    f2c[i][2]  = sigma(x_Bj, q2charm[i],0.0,par_temp);
	}

	printf("%.3e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n",x_Bj,
                f2c[0][0]/f2a[0][0],f2c[0][1]/f2a[0][1],f2c[0][2]/f2a[0][2],
                f2c[1][0]/f2a[1][0],f2c[1][1]/f2a[1][1],f2c[1][2]/f2a[1][2],
                f2c[2][0]/f2a[2][0],f2c[2][1]/f2a[2][1],f2c[2][2]/f2a[2][2],
                f2c[3][0]/f2a[3][0],f2c[3][1]/f2a[3][1],f2c[3][2]/f2a[3][2]);
	fprintf(fout,
                "%.3e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n",x_Bj,
                f2c[0][0]/f2a[0][0],f2c[0][1]/f2a[0][1],f2c[0][2]/f2a[0][2],
                f2c[1][0]/f2a[1][0],f2c[1][1]/f2a[1][1],f2c[1][2]/f2a[1][2],
                f2c[2][0]/f2a[2][0],f2c[2][1]/f2a[2][1],f2c[2][2]/f2a[2][2],
                f2c[3][0]/f2a[3][0],f2c[3][1]/f2a[3][1],f2c[3][2]/f2a[3][2]);

    }
    fclose(fout);
}



/*******************************************************************************
* The function plotting dipol cross section 
*******************************************************************************/
void graph_f_2_beauty (int model_graph) {

    FILE *fout; 
     int i,k;

    double q2charm[6] = {3.5,12,25,60,200,650};
    double x_Bj, f2b,f2b_mod;
    double par_temp[5];

    double npoints; /* Number of points used to plot */
    ////////////////double r = 0.0;

    model = model_graph;

    switch (model) {
    case 0:
	sigma_0 = 29.12; 
        lambda  = 0.277; 
        x_0     = 0.41e-04; 
/*
	sigma_0 = par_start[0]; 
        lambda  = par_start[1]; 
        x_0     = par_start[2]; 
*/
	par_temp[0] = sigma_0; 
        par_temp[1] = lambda; 
        par_temp[2] = x_0; 
        printf("GBW\n");
	fout = fopen("f2_beauty_gbw.dat","w");    
    break;
    case 1:
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

        printf("BGK\n");
	fout = fopen("f2_beauty_bgk.dat","w");    
	chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);
    break;
    }

    npoints = 10;
    flavour = 3;

	//fout = fopen("graph_f_2_beauty.dat","w");    

	for(i=0;i<2;i++) {
            for(k=0;k<npoints;k++) {
                //x_Bj = 0.01; /* x = 10^-1...10^-6 */
                x_Bj = pow(10,-(0.5*k+2)); /* x = 10^-1...10^-6 */
		f2b     = sigma(x_Bj, q2charm[i],0.0,par_temp);
		f2b_mod = pow((1-x_Bj),20)*f2b;
	        fprintf(fout,"%e %e %e  %e\n", 
                      q2charm[i], x_Bj, f2b,f2b_mod);
            }
	    fprintf(fout,"\n"); 
        }

	for(i=2;i<6;i++) {
            for(k=0;k<npoints;k++) {
                //x_Bj = 0.01; /* x = 10^-1...10^-6 */
                x_Bj = pow(10,-(0.5*k+1)); /* x = 10^-1...10^-6 */
		f2b     = sigma(x_Bj, q2charm[i],0.0,par_temp);
		f2b_mod = pow((1-x_Bj),20)*f2b;
	        fprintf(fout,"%e %e %e  %e\n", 
                      q2charm[i], x_Bj, f2b,f2b_mod);
            }
	    fprintf(fout,"\n"); 
        }

	fclose(fout);

    //update_graph("f_2_charm.plt");
}

/*******************************************************************************
* The function plotting F_2_charm
*******************************************************************************/
void graph_dipolcs (void) {

    FILE *fout; 
     int i,k;

    double npoints; /* Number of points used to plot */
    double r, rmin = 0.01;
    double delta;
    double r_list[6] = {0.1,0.5,1,5,10,20};
    double dcsr[5],dcsx[5];

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

    npoints = 50;
/*
    xmod = 0.001;
    printf("%f %f %f \n",sigma_0,lambda,x_0);
    printf("%f\n",sigma_gbw(1.0));
*/
    delta = (log10(30)-log10(0.01))/npoints; /* Logarithmic x-scale */

	fout = fopen("dipol_csr.dat","w");    
	for(i=0;i<=npoints;i++) {
            r = rmin*pow(10,delta*i);
            for(k=0;k<5;k++) {
                xmod = pow(10,-(k+2)); /* x = 10^-2...10^-6 */
		switch (model) {
		    case 0:
		    dcsr[k] = sigma_gbw(r);
		    break;
		    case 1:
		    dcsr[k] = sigma_bgk(r);
		    break;
                } 

            }

	    fprintf(fout,"%e %e %e %e %e %e\n",
                    r, dcsr[0], dcsr[1], dcsr[2], dcsr[3], dcsr[4]);
        }
	fclose(fout);

    delta = (log10(1.0e+00)-log10(1.0e-07))/npoints; /* Logarithmic x-scale */
	fout = fopen("dipol_csx.dat","w");    
	for(i=0;i<=npoints;i++) {
            xmod = 1.0e-07*pow(10,delta*i);
            for(k=0;k<6;k++) {
                r = r_list[k];
                switch (model) {
		    case 0:
		    dcsx[k] = log(sigma_0 - sigma_gbw(r));
		    break;
		    case 1:
		    dcsx[k] =  sigma_bgk(r);
		    //dcsx[k] = log(sigma_0 - sigma_bgk(r));
		    break;
                } 
                
            }

	    fprintf(fout,"%e %e %e %e %e %e %e\n",
                    xmod, dcsx[0], dcsx[1], dcsx[2], dcsx[3], dcsx[4], dcsx[5]);
        }
	fclose(fout);

    update_graph("dipol_cs.plt");
}
/*******************************************************************************
* The function printing info with starting values
*******************************************************************************/
void printinfo (char *out_name) {

    FILE *out_file;
    char *fl;

    if(strcmp(out_name,"stdout")==0)
        out_file = stdout;
    else 
	out_file = fopen(out_name,"w");

    if (flavour==0) fl = "l";
    else if (flavour==1) {
        if (fl_beauty==0) fl = "lc";
        else if (fl_beauty==1) fl = "lcb";
    }

    printf("\n");
    fprintf(out_file,"Started...\n");
    fprintf(out_file,"-------------------------------------------------------");
    fprintf(out_file,"-------------------------\n");
    fprintf(out_file,
           " FLAVOUR:   %s       MX:  %d     eps_simps2d:   %1.1e",   
            fl,MX, eps_simps2d);
    fprintf(out_file,
           "    n_f:     %d\n",n_f);
    fprintf(out_file," ACTION:    %s     MQ:  %d     Rmax_simps2d:  %1.1e",\
            (action==0)?"fit ":"list", MQ, Rmax_simps2d);
    fprintf(out_file,
           "    Q^2_0:   %1.1f\n",Q0);
    //char modelname[10];
    //if(model==0){modelname="GBW"; }else if(model==1){modelname= "BGK";}else if(model==2){modelname="GBS";}
    fprintf(out_file," MODEL:     %s                      m_l[GeV]:  %.2f",
		  (model==0)?"GBW":((model==1)?"BGK":((model==2)? ( (sudflag==0)?"GBS, no S":((sudflag==1)?"GBW, Spert":((sudflag==2)?"GBW, S": "N/A") )):"N/A")),
		    sqrt(m_fsq)
	);
            //(model==0)?"gbw":"bgk",sqrt(m_fsq));
    fprintf(out_file,
           "     q_down:   %1.1e\n",q_down);
    fprintf(out_file," DATAFORM:  %s                      m_c[GeV]:  %.3f ",
            (dataform==0)?"red_cs":"F_2",sqrt(m_ch));
    fprintf(out_file,
           "      q_up:   %1.1e\n",q_up);
    fprintf(out_file," labda_qcd[MeV]: %.0f                 m_b[GeV]:  %.3f\n",
            sqrt(Lambda2)*1000, sqrt(m_b));
    fprintf(out_file,"-------------------------------------------------------");
    fprintf(out_file,"-------------------------\n");
    fprintf(out_file,"\n");
    
    //if((out_name != stdout))
    if(strcmp(out_name,"stdout")!=0)
	fclose(out_file);
    
}

/*******************************************************************************
* 
*******************************************************************************/
void save_result (void) {

    FILE *result_file;

    result_file = fopen("fit_result.txt","w");

    fprintf(result_file,"%e %e %e %e %e\n",sigma_0,A_g,lambda_g,C,mu02);
    fclose(result_file);


}

/*******************************************************************************
* The FCN function
*******************************************************************************/
void fcn (int npar, double grad[], double *fcnval,
          double *par,int iflag, void (*Dummy)()) {
          
    double chipar;  
    int    nbins;  /* Number of points used to calculate chi^2 */       
    int        i; 


    //printf("HERE\n");
    /* Initial value of chi^2 equal to '0' */
    *fcnval=0.0;  

    /* Set number of bins */
    switch(flavour) {   
    case 0: /* Light  only */
    case 1: /* Light + charm */
        nbins = nf2data;
        //nbins = NDATA;
        break;
    case 2: /* Charm only */
        break;
    }

    switch(model) {
    case 0:
	sigma_0  = par[0]; 
	lambda   = par[1];
	x_0      = par[2];
    break;
    case 1:
	sigma_0  = par[0]; 
	A_g      = par[1];
	lambda_g = par[2];
	C        = par[3]; 
	mu02     = par[4];
 
	chebft(xmin,xmax,Qmin,Qmax,MX,MQ,coef, &xgpdf);
    break;
    case 2:
	sigma_0  = par[0]; 
	lambda   = par[1];
	x_0      = par[2];
        C        = par[3]; 
        //mu02     = par[4];
	r_max 	=par[4];
//	bmax	= par[4];
	g1	 = par[5];
	g2 	 = par[6];
  
        //chebft3(xmin,xmax,Qmin,Qmax,rmin,rmax,NX,NQ,NR,coef3,&logS_gbs);
        //printf("chebft3 done...\n", nbins);
	//exit(1);

    break;
    }
     //fprintf(stdout,"Initial: %e %e %e %e %e %e %.4f\n",par[0],par[1],par[2],par[3], par[4],par[5] ,*fcnval/nf2data);
    //printf("%d\n", nbins);

    /* Calcualte chi^2 */
    for(i=0;i<nbins;i++) {
       sigmafcn[i]=sigma(xmoddata[i],q2data[i],ydata[i],par);
       //printf("%e %e %f\n",xmoddata[i],q2data[i],sigmafcn[i]);
       chipar = (csdata[i]-sigmafcn[i])/erdata[i];
       //printf("%f\n",sigmafcn[i]);
       //printf("%f\n",chipar);
       //printf("%f\n",csdata[i]);
       //printf("%e %e %e %e %e %e \n", 
               //par[0], par[1], par[2], par[3], par[4], chipar);
        *fcnval+=chipar*chipar;
        //printf("%3d  %5.1f  %e  %.8e  %.8e  %.8e\n",
        //        i+1,q2data[i], xmoddata[i],erdata[i],csdata[i],sigmafcn[i]);
        //printf("%e\n",*fcnval);
    }
    /* Print current parameters values and chi^2 */

   printf("%.3e %.3e %.3e ",par[0],par[1],par[2]);
   if(model>=1)printf("%.3e %.3e ",par[3], par[4]);
   if(sudflag>=2)printf("%.3e %.3e ",par[5],par[6]);
   //if(sudflag>=1)printf("%.3e ",par[6]);
   printf("%.4f /%.4f = ",(*fcnval) ,(float) nf2data);
   printf("%.4f",((*fcnval)/((float)nf2data)) );
   printf("\n");
    //graphdata();
    //exit(1);
}

