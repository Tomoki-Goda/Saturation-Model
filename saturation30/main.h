
/* The number of proper experimental data for total reduced cross section; 
   the number of all data is 331 */
#define NDATA 500
//#define  NDATA 318   // 50
//#define  NDATA 222   //10

/* Global arrays with the values of experimental data */
double   xbdata[NDATA];  /* Bjorken variable xB data*/
double    ydata[NDATA];
double    wdata[NDATA];
double   q2data[NDATA];
double   csdata[NDATA];
double   erdata[NDATA];
double sigmafcn[NDATA];
double    sfcn0[NDATA];
double    sfcn1[NDATA];

/* Global array with the values of modified Bjorken variable 
 *         +----------------------------------------------+
 *         |                                              |
 *         |  xmoddata[i] = xbdata[i](1+4m_fsq/q2data[i]) |
 *         |                                              |
 *         +----------------------------------------------+ 
 */
double xmoddata[NDATA];

/* Limits for Q^2 and Bjorken x used to fit */
double q_down = 0.0;
double q_up   = 10.0;
double x_up   = 0.01;

int nf2data;
/* Rescale factor for H1 data */
double rescale = 1.05;


double q2graph[25]={1.5,2.0,2.5,3.5,5.0,6.5,8.5,10,12,15,18,22,25,35,45,60,70,
                    90,100,120,150,200,250,300,400};

/* Global variables */
//double arg1,arg2;
double   Q2;
double   xmod1; 
double   xmod;   /* Modified Bjorken variable 
                  *   +---------------------------+ 
                  *   |   xmod = xB(1+4m_fsq/Q2)  | 
                  *   +---------------------------+ 
                  */
double   y;

/*Model parameters */
double sigma_0;
/* GBW */
double  lambda;
double     x_0;
/* BGK */
double  A_g;
double  lambda_g;
double  C;
double  mu02;

/*GBS-G*/
double g1, g2;


/* Initial condition */
//double       beta = 9.6;
double       beta = 6.6;

/* GBW Starting parameter values, errors: sigma_0,lambda,x_0, C, mu2,g1 */
double par_start[6]={ 23.0 ,0.29 ,3.0e-4 ,1.26 ,5.0, 0.2};
double par_error[6] = { 1.0,  0.05,  0.1e-04 ,0.01,0.01, 0.01 };

double   par_min[6] = { 0.0,  0.00,  0.0,0.01, 0.0 ,  0.0};
double   par_max[6] = {80.0,  1.00,  1.0,20.0, 20.0,  1.0 };

/* BGK Starting parameter values, errors: sigma_0, A_g, lambda_g, C, mu02 */

double bgk_parst[5]  = {22.40,1.35,0.079,0.38, 1.73}; /**/
double bgk_parer[5]  = { 1.00,  0.10,  0.05, 0.05, 0.40};
double bgk_parmin[5] = { 0.00,  0.00,  0.01, 0.01, 0.50};
double bgk_parmax[5] = {30.00,  5.00,  0.20, 1.00, 5.00};

/* Precision parameters */
int    simpsN       = 500;  /* was 400 */
double sigmaEPS     = 1.0e-06;
double eps_simps2d  = 1.0e-02; /* 1.0e-05  -  does not improve much 
                                  1.0e-03  -  still very accurate !  */
double STRategy     = 0.0;
double Rmax;
double amin  = 1.0e-06;
double Rmax_dadmul  = 3.0e+01;
//double Rmax_dadmul  = 5.0e+03;
double Rmax_simps2d = 3.0e+01;
//double Rmax_simps2d = 1.0e+04;


/* Constants */
double         pi = 3.141592653589793;
double         PI = 3.141592653589793;

double eulergamma = 0.577215664901532;

double      m_fsq = 0.0196;  /* Common light quark mass in GeV^2 */
double       m_s  = 0.0196;  /* Strange quark mass in GeV^2 */
double       m_ch = 1.96;    /* Cham quark mass in GeV^2 */
double       m_b  = 21.16;   /* Bottom quark mass in GeV^2 */

double       norm = 0.08282813;     /* \Psi^2 normalization */
double       hc2  = 0.389379;
int           n_f = 3;         /* Number of avtive flavours */
double        n_0 = 0.5;       /* Maximal singluraity of integrand */
double         Q0 = 1.0;    
double    Lambda2 = 0.09;    /* QCD Lambda^2 = (300)^2 GeV^2 */
double	bmax=0.5;


/* CERNLIB functions */
extern  double dadmul_();  /* Adaptive quadrature integral */
extern  double dbesk0_();  /* Modified Bessel Function K_0 */
extern  double dbesk1_();  /* Modified Bessel Function K_1 */
extern  double dbesj0_();  /* Bessel Function J_0 */
extern  double dbesj1_();  /* Bessel Function J_1 */
//extern  double dgquad_(); 
double simps2d(double *xmi, double *xma,
               double precision, double (*func)(int*,double*));

double dbeta;
double W2;
double B_D = 6.0;
double xp;
double negglu;
