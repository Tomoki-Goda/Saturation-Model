/*  
 *  Claculating of xgpdf(x,Q^2) by DGLAP evolution of initial condition 
 *  xgpdf(x,Q^2_0)
 * 
 * Originally written by S. Sapeta 
 * Modified by T. Goda
 */

//#include "main.h"
#include "../complex.h"


////////////// use values from the control & constant ///////////////////
static double n_f=NF;
static double Lambda2=LQCD2;
static double pi=PI;
static double eulergamma = GAMMA_E;
////from main.h////////
static double simpsN=500;
/* Initial condition */
//double       beta = 9.6;
double       beta = 6.6;
double        n_0 = 0.5;       /* Maximal singluraity of integrand */

////from control.h/////////////
int  gluon_int    = 1;
//////////////////       fit parametr          /////////////////////////
static double A_g, lambda_g;
////////////////////////////////////////////////////////////////////////

///////////////////function to control fitparameter from outside the file////////////////
void set_xg_parameter(double ag,double lg){
	A_g=ag;
	lambda_g=lg;
}
//////////////////////////////////////////////////////////////////////////////////////////


////////////// Old Part, no major change below/////////////////////////////////////////

/* CERNLIB functions*/
extern doublecomplex wgamma_();
extern doublecomplex wpsipg_();
extern        double dgammf_();
extern        double dsimps_();
extern        double simps_();

double Yg, tg;

/////////////////////////////////////////////////////////////////////////////////////
//////////////////   GLUON DISTRIBUTION xgpdf(x,Q^2)  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


/* initial condition */
double xgpdf_init(double x) {
	return A_g*pow(x,-lambda_g)*pow((1-x),(beta-1));
}

/* gammatilde function */
doublecomplex gammatilde(doublecomplex n) {

    doublecomplex n1,n2,n3,l1,l2,t1,t2,t3,cx,value;
           double m1,m2,rl;
              int k=0;

    n1 = Cadd(n,Complex(1,0));
    n2 = Cadd(n,Complex(2,0));
    n3 = Cadd(n,Complex(3,0));

    l1 = Cmul(Conjg(n),Conjg(n1));
    m1 = Cabs(n)*Cabs(n)*Cabs(n1)*Cabs(n1);
    t1 = RCmul(1/m1,l1);

    l2 = Cmul(Conjg(n2),Conjg(n3));
    m2 = Cabs(n2)*Cabs(n2)*Cabs(n3)*Cabs(n3);
    t2 = RCmul(1/m2,l2);

    t3 = wpsipg_(&n2,&k);

    cx = Csub(Cadd(t1,t2),t3);
    rl = 11.0/2.0-n_f/3.0-6.0*eulergamma;

    value = Cadd(RCmul(6,cx),Complex(rl,0));
    return value;
}


/* integrand of gluons pdf */
double xgpdf_integrand (double *y) {
//doublecomplex xgpdf_integrand(double y, double Y, double t) {

    doublecomplex n0,n1,n2,g1,g2,gt,ex,l;
    //doublecomplex n0,n1,n2,g1,g2,gt,ex,l,value;
           double val;
           double m;

    //n0 = Complex(n_0,y);
    //n1 = Complex(-lambda_g+n_0,y);
    //n2 = Complex(-lambda_g+beta+n_0,y);
    n0 = Complex(n_0,y[0]);
    n1 = Complex(-lambda_g+n_0,y[0]);
    n2 = Complex(-lambda_g+beta+n_0,y[0]);

    g1 = wgamma_(&n1);
    g2 = wgamma_(&n2);

    //gt = RCmul(t,gammatilde(n0));
    gt = RCmul(tg,gammatilde(n0));

    //ex = Cexp(Cadd(Complex(0,y*Y),gt));
    ex = Cexp(Cadd(Complex(0,y[0]*Yg),gt));

    l = Cmul(g1,Conjg(g2));
    m = Cabs(g2)*Cabs(g2);

    //value = Cmul(RCmul(1/m,l),ex);
    val = Cmul(RCmul(1/m,l),ex).r;
    return val;
    //return value;
}

doublecomplex xgpdf_integrand2(double y, double Y, double t) {

    doublecomplex n0,n1,n2,g1,g2,gt,ex,l,value;
           double m;

    n0 = Complex(n_0,y);
    n1 = Complex(-lambda_g+n_0,y);
    n2 = Complex(-lambda_g+beta+n_0,y);

    g1 = wgamma_(&n1);
    g2 = wgamma_(&n2);

    gt = RCmul(t,gammatilde(n0));

    ex = Cexp(Cadd(Complex(0,y*Y),gt));

    l = Cmul(g1,Conjg(g2));
    m = Cabs(g2)*Cabs(g2);

    value = Cmul(RCmul(1/m,l),ex);
    return value;
}
/*
double integrand(double *y) {

    return xgpdf_integrand(y[0], Yg, tg);
}
*/

void xgpdf_av(void) {

    double av;
    double a1,ab1;
    
//    lambda_g = -0.0908;
//    A_g = 2.368;

    a1 = -lambda_g +1;
    ab1 = -lambda_g + 1 + beta;

    av = A_g*dgammf_(&a1)*dgammf_(&beta)/dgammf_(&ab1);
    //value = 0.2;

   printf("<xgpdf(x,Q^2_0)> = %f\n",av);
}

/*******************************************************************************
* Gluons pdf  function
*******************************************************************************/
double xgpdf(double x, double QQ) {
          
     /****************************************************\
     *  Integration may be performed using two methods:  *
     *                                                   *
     *   -   Golec Simpson method  (gluon_int = 0)       *
     *   -   CERN Simpson method   (gluon_int = 1)       *
     *                                                   *
     \****************************************************/

    double a = 0.0, b;
    double normalization;
    double value;
    double bprim; 

    /* Golec Simpson variables */
    double c = 150.0;
    double prec = 1.0e-08;
    double dum1,dum2,dum3;

    /* CERNLIB Simpson variables */
    int    N; // must be even !!!
	   N = simpsN; /* Precision of integration control by global
			  variable declared in main.h. N is the number
			  of points in which the integrand function
			  is calculated.                              */
			
    int    i;
    double fr[N+1];

    bprim = 33.0/6.0-n_f/3.0;
    Yg = log(1/x);
    tg = (1/bprim)*log(log(QQ/Lambda2)/log(Q0/Lambda2));

    normalization = A_g*exp(n_0*Yg)*dgammf_(&beta)/pi;
    
    switch (gluon_int) {
        /* Golec Simpson method of 1d integration */
        case 0:
        simps_(&a,&c,&prec,&prec,&prec,
               &xgpdf_integrand,&dum1,&value,&dum2,&dum3);   
        break;

        /* CERNLIB Simpson method of 1d integration */
        case 1:
        b = 180/(log10(QQ)+3);  /* Upper limit of integration */ /* was 30 */

        for(i=0;i<=N;i++) {
	         fr[i]=xgpdf_integrand2(i*(b-a)/N,Yg,tg).r;

        }
	//printf("simpsN %d\n",N);
        value = dsimps_(fr,&a,&b,&N);
        break;
    }
    return  normalization*value;
}



