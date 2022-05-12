#include  <stdio.h>
#include  <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "control.h"
#include "test.h"
#include "cfortran.h"
#include "minuit.h"

double ftest(double r) {
  return r+r*r;
}

//double logS_gbs(double r) {
//  double val = S_gbs(r);
//  //printf("%f %e\n", r, val);
//  
//  if(val<1e-50) val = -log(1e-50);
//  else val = -log(val);
//  return val;
//}

// GBWS update
//double sigma_0 = 23.58;
//double lambda = 0.27;
//double x_0 = 2.24;

// GBS update
double sigma_0 = 21.39;
double lambda = 0.26;
double x_0 = 2.23e-4;
double C = 2.19;

int main(int argc, char *argv[]) {

  //chebft(xmin,xmax,rmin,rmax,NX,NR,coef2,&logS_gbs);

  double N = 10;
  double r1min = 1e-3;
  xmod = 0.00001;
  Q2= 2.0;
  for (int j=0; j<N; j++) {
    double r1 = pow(rmax/r1min,j/N)*r1min;
    //printf("%e %e\n", r1, sigma_gbs(r1));
    printf("%e %e %e\n", r1, sigma_gbs(r1), sigma_gbs_cheb(r1));
    //printf("%e %e %e\n", r1, sigma_gbs(r1), sigma_gbs_cheb(r1));
    //printf("%e %e %e\n", r1, sigma_gbw(r1), sigma_gbs(r1));
  }

  //printf("%e\n", S_gbs(xmod,1e-30));
  //printf("%e %e\n", 1e-70, S_gbs(xmod,1e-70));
  //printf("%e %e\n", 1e-60, S_gbs(xmod,1e-60));

  //for (int j=0; j<N; j++) {
  //  double r1 = pow(rmax/0.001,j/N)*r1min;
  //  //printf("%e %e\n", r1, exp(-sudakov(r1)));
  //  printf("%e %e\n", r1, S_gbs(xmod,r1));
  //  //printf("%e %e  %e\n", r1, sigma_gbw(r1), sigma_0*(1-S_gbs(xmod,r1)));
  //}

   //chebft(xmin,xmax,rmin,rmax,NX,NR,coef2,&logS_gbs);

   //double x1 = 1.702920e-06; 
   //double r1 = 1.000000e-03;

   //xmod = x1;
   //printf("%e %e %e %e  %e\n", x1, r1, sigma_gbw(r1), sigma_gbs(r1),
   //       sigma_0*(1-exp(-chebev(xmin,xmax,rmin,rmax,NX,NR,coef2,x1,r1))));
                               //S_gbs_cheb(xmod, r1));
                               //sigma_0*(1-S_gbs_cheb(xmod, r1)));
   //printf("%e %e %e %e\n", x1, r1, 1-exp(-logS_gbs(x1, r1)), 
   //                       1-exp(-chebev(xmin,xmax,rmin,rmax,NX,NR,coef2,x1,r1)));
  //////printf("%e %e\n", rmin, rmax); 
  //double N = 20;
  //for (int i=0; i<N; i++) {
  //  double x = pow(xmax/xmin,i/N)*xmin;
  //  //double x = (xmax-xmin)*i/N+xmin;
  //  for (int j=0; j<N; j++) {
  //    //double r = (rmax-rmin)*j/N+rmin;
  //    double r = pow(rmax/0.001,j/N)*0.001;
  //    printf("%e %e %e %e\n", x, r, logS_gbs(x, r), 
  //                         chebev(xmin,xmax,rmin,rmax,NX,NR,coef2,x,r));
  //  }
  //}


  //chebft1(rmin,rmax,MR,coef1, &logS_gbs);
  ////chebft1(rmin,rmax,MR,coef1, &S_gbs);

  //for (int i=0; i<1000; i++) {
  //  double r = (10-0.001)*i/100+0.001;
  //  if (r>rmax) break;
  //  //printf("%f %e %e\n", r, logS_gbs(r), chebev1(rmin,rmax,MR,coef1,r));
  //  //printf("%f %e\n", r, 1-S_gbs(r));
  //  
  //  //printf("%e %e %e\n", r, sigma_gbs(r), sigma_gbw(r));

  //  printf("%e %e %e\n", r, logS_gbs(r), chebev1(rmin,rmax,MR,coef1,r));
  //  //printf("%f %e %e\n", r, 1-exp(-logS_gbs(r)), 
  //  //       1-exp(-chebev1(rmin,rmax,MR,coef1,r)));

  //  //printf("%f %e %e\n", r, logS_gbs(r), chebev1(rmin,rmax,MR,coef1,r));
  //  //printf("%f %e %e\n", r, log(S_gbs(r)), chebev1(rmin,rmax,MR,coef1,r));
  //  //printf("%f %e %e\n", r, -log(S_gbs(r)), chebev1(rmin,rmax,MR,coef1,r));
  //}

  //r = 0.3;
  //printf("%f %f\n", ftest(r), chebev1(rmin,rmax,MR,coef1,r));
  //r = 1.3;
  //printf("%f %f\n", ftest(r), chebev1(rmin,rmax,MR,coef1,r));
  //r = 0.001;
  //printf("%f %f\n", ftest(r), chebev1(rmin,rmax,MR,coef1,r));
  //r = 1300;
  //printf("%f %f\n", ftest(r), chebev1(rmin,rmax,MR,coef1,r));

  //double xi, yi, x[10], y[10];
  //for (int i = 0; i < 10; i++) {
  //    x[i] = i;
  //    y[i] = ftest(i);
  //}
  //gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  //gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, 10);
  //gsl_spline_init (spline, x, y, 10);

  //for (xi = x[0]; xi < x[9]; xi += 0.1) {
  //  yi = gsl_spline_eval (spline, xi, acc);
  //  printf ("%g %g %g\n", xi, yi, ftest(xi));
  //}

  return 0;
}
