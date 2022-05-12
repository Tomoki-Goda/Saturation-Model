#include  <stdio.h>
#include  <math.h>
#include "control.h"
#include "cfortran.h"
#include "functions.h"

//double func(double x, double y, double z) { return x+y/z; }

int main(int argc, char *argv[]) {


  FILE *parfile;
  parfile = fopen(argv[1],"r");
  fscanf(parfile, "%lf %lf %lf %lf", &sigma_0,&lambda,&x_0,&C);

  //printf("%d\n",strcmp(argv[1],"par1.dat"));
  //printf("#%e %e %e\n", sigma_0, lambda, x_0);

  xmod = 1e-3;
  Q2 = 100.0;

  double N = 100;
  double r1min = 1e-3;


  if (strcmp(argv[1],"par0.dat")==0) {
    for (int j=0; j<N; j++) {
      double r1 = pow(rmax/r1min,j/N)*r1min;
      printf("%e %e\n", r1, sigma_gbw(r1));
    }
  } else {
    //chebft3(xmin,xmax,Qmin,Qmax,rmin,rmax,NX,NQ,NR,coef3,&logS_gbs);
    for (int j=0; j<N; j++) {
      double r1 = pow(rmax/r1min,j/N)*r1min;
      printf("%e %e\n", r1, sigma_gbs_exact(r1));
    }
  }

  return 0;
}
