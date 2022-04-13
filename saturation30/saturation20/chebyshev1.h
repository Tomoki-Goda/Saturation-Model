//
//     Author: Sebastian Sapeta
//       Date: August 2005
//

#define  MR 20   /* The order of polynomial in x variable */

double coef1[MR];     /* The coefficients vector */
double rmin =  1.0e-06;  /* x lower limit */
double rmax =  3.0e+01;  /* x upper limit */
// Too large rmax leads to worse fit of the low-r part of the function and this
// has an effect on the integral. On the contrary taking chebeshev at r>rmax
// in the integral does not affect the result as this part of dipole cross
// section is in any case suppreseed.


/*******************************************************************************
* Chebyshev - coefficients calculation
*******************************************************************************/
extern void chebft1 (double rmin, double rmax, int n,
             double c[MR] , double (*func)(double)) {
    int    k,i;  /* Counters */
    double fac = 2.0/n;
    double sum;
    double f[MR]; /* Function values array */ 

    double yr; /* Polynomial's zeros in y(x)*/
    double zerox; /* Polynomial's zeros in x */

    /* Calculation of function values at zeros */
    for (k=0;k<n;k++) {
            /* Calculating zeros in y(x) and y(Q^2) */
            yr=cos(pi*(k+0.5)/n); 

            /* Calculating zeros in x and Q^2 */
            zerox = 0.5*(rmax-rmin)*yr+0.5*(rmax+rmin);
            //zerox = rmax*pow((rmax/rmin),(0.5*(yr-1)));

            /* Function values at zeros */
            f[k]=func(zerox);  
	    //printf("%e %e  %e\n", zerox, f[k]);
    }


    /* Calculation of Chebeshev coefficients  */
    for (i=0;i<n;i++) {
            sum = 0.0;
            for(k=0;k<n;k++) { sum +=f[k]*cos(pi*i*(k+0.5)/n); }
        c[i]=fac*sum; 
   }
}


/*******************************************************************************
* Chebyshev evaluation - returns the value of approximated function
*******************************************************************************/
double chebev1 (double rmin, double rmax, int n,  double c[MR] , double r) {

    double yrcheb[MR];

    /* Transformation coefficients */
    double Ar = 2.0/(rmax-rmin);
    double Br = -(rmax+rmin)/(rmax-rmin);
    //double Ar =  2/log(rmax/rmin);
    //double Br =  1-2/log(rmax/rmin)*log(rmax);
    //double Ar = 2.0/log(rmax/rmin);
    //double Br =  1-2.0/log(rmax/rmin)*log(rmax);
    //double Br =  1;
     
    double yr;
    double sumr = 0.0;/* Partial summs */
    int i,k; /* Counters */

    /* Logarithmic transformation*/
    //yr = r;
    yr = Ar*r+Br;
    //yr = (Ar*log(r)+Br);

    //printf("%f\n",yr);

    yrcheb[0] = 1;
    yrcheb[1] = yr;
    
    for (k=2;k<MR;k++) yrcheb[k]=2*yr*yrcheb[k-1]-yrcheb[k-2];   

    for (i=1;i<n;i++) sumr += c[i]*yrcheb[i];
	 
    return 0.5*c[0]+sumr;

}
