#include  <stdio.h>
#include  <math.h>
#include "control.h"
#include "cfortran.h"
#include "functions.h"

//double func(double x, double y, double z) { return z; }
double func(double x, double y, double z) { return exp(-z); }

int main(int argc, char *argv[]) {

  // GBS update
  sigma_0 = 1.0;
  //sigma_0 = 21.39;
  lambda = 0.26;
  x_0 = 2.23e-4;
  C = 0.5;

  double N = 100;
  double r1min = rmin;
  //double r1min = 1e-3;
  double q2min = 5e-1;
  xmod = 1e-4;
  Q2 = 1.0;

  //chebft3(xmin,xmax,Qmin,Qmax,rmin,rmax,NX,NQ,NR,coef3,&logS_gbs);
  //chebft3(xmin,xmax,Qmin,Qmax,rmin,rmax,NX,NQ,NR,coef3,&func);


  for (int j=0; j<N; j++) {
    double r1 = pow(rmax/r1min,j/N)*r1min;
    //printf("%e %e %e\n", r1, func(r1, r1, r1),
    //       chebev3(xmin,xmax,Qmin,Qmax,rmin,rmax,NX,NQ,NR,coef3,r1,r1,r1));
    //printf("%e %e %e\n", r1, func(xmod, Q2, r1),
    //       chebev3(xmin,xmax,Qmin,Qmax,rmin,rmax,NX,NQ,NR,coef3,xmod,Q2,r1));
    //printf("%e %e %e\n", r1, S_gbs(xmod, Q2, r1), S_gbs_cheb(xmod, Q2, r1));
    //printf("%e %e %e\n", r1, logS_gbs(xmod, Q2, r1),
    //       chebev3(xmin,xmax,Qmin,Qmax,rmin,rmax,NX,NQ,NR,coef3,xmod,Q2,r1));
    //printf("%e %e %e\n", r1, sigma_gbs_exact(r1), sigma_gbs(r1));
    printf("%e %e %e\n", r1, sigma_gbs_exact(r1), sigma_gbw(r1));
  }

  //double r1 = 0.3;
  //for (int j=0; j<N; j++) {
  //  double q2 = pow(Qmax/q2min,j/N)*q2min;
  //  printf("%e %e %e\n", q2, logS_gbs(xmod, q2, r1),
  //         chebev3(xmin,xmax,Qmin,Qmax,rmin,rmax,NX,NQ,NR,coef3,xmod,q2,r1));
  //}


  return 0;
}
