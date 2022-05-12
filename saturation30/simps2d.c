/*  
 *    Author: Sebastian Sapeta
 *      Date: October 2005
 *
 */

#include     <stdio.h>
#include     <math.h>

/* 1-d integration by Simpson method */
extern double simps_();
extern double dgquad_();


static double arg1;     /* Second argument of function */
//static double arg2;     /* Second argument of function */
static double eps;      /* Precision of integration */
static double argmin[2];  /* Lower limits */
static double argmax[2];  /* Uper limits */
static double dum1_f1, dum2_f1, dum3_f1;  /* Dummy variables */
static double result;   /* Result of integration */

/* Function of two variables */
//static double (*f)(double,double);
static double (*f)(int*,double*);


/* First integration function */
//double f1 (double x) {
double f1 (double *x) {
   //double arg1 = x[0];
   //double arg2 = x;
   double arg2 = x[0];
   int d = 2;
   double xtemp[2] = {arg1,arg2};
   //return f(arg1,arg2);

   //return 0.0;
   return f(&d,xtemp);
}


/* First integration function - over arg1 i.e. r */
double f2 (double *x) {

    double value;
    double dum1_f2,dum2_f2,dum3_f2; /* Dummy variables */
    //int N = 64;
    int N = 96;


    arg1 = x[0];   /* First integration function - over arg1 i.e. z */
    //arg2 = x[0];



    value = dgquad_(&f1,&argmin[0],&argmax[0],&N);

    //simps_(&argmin[0],&argmax[0],
    //       &eps,&eps,&eps,&f1,&dum1_f2,&value,&dum2_f2,&dum3_f2);

    return value;
}



/*******************************************************************************
* 2-d integral of function by Simpson method 
*******************************************************************************/
double simps2d(double *xmi, double *xma, 
               double precision, double (*func)(int*,double*)) {

    /* Precision */
    eps = precision;

    /* Limits */
    argmin[1] = xmi[0];
    argmin[0] = xmi[1];
    argmax[1] = xma[0];
    argmax[0] = xma[1];
/*
    argmin[0] = xmi[0];
    argmin[1] = xmi[1];
    argmax[0] = xma[0];
    argmax[1] = xma[1];
*/
    /* Integrated function */
    f = func; 

    /* Second integration function - over arg1 i.e. z */
    simps_(&argmin[1],&argmax[1],
	   &eps,&eps,&eps,&f2,&dum1_f1,&result,&dum2_f1,&dum3_f1);

    return result;
}

