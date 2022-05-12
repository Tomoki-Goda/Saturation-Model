//
//     Author: Sebastian Sapeta
//       Date: August 2005
//

#define  NX 10   /* The order of polynomial in x variable */
#define  NQ 10   /* The order of polynomial in x variable */
#define  NR 20  /* The order of polynomial in x variable */
//double rmin =  2.0e-00;  /* r lower limit */
//double rmax =  2.0e+01;  /* r upper limit */
double rmin =  1.0e-06;  /* r lower limit */
double rmax =  3.0e+01;  /* r upper limit */
double coef3[NX][NQ][NR];     /* The coefficients vector */

//double  pi = 3.141592653589793;

/*******************************************************************************
* Chebyshev - coefficients calculation
*******************************************************************************/
extern void chebft3 (double xmi, double xma, 
             double Qmi, double Qma, 
             double rmi, double rma, 
             int n, int m, int o,
             double *c , double (*func)(double,double, double)) {
    int    k,l,q,i,j,r;  /* Counters */
    double fac = 8.0/(n*m*o);
    double sum;
    double f[n][m][o]; /* Function values array */ 

    double yx; /* Polynomial's zeros in y(x)*/
    double yQ; /* Polynomial's zeros in y(Q^2) */
    double yr; /* Polynomial's zeros in y(Q^2) */

    double zerox; /* Polynomial's zeros in x */
    double zeroQ; /* Polynomial's zeros in Q^2 */
    double zeror; /* Polynomial's zeros in Q^2 */



    /* Calculation of function values at zeros */
    for (k=0;k<n;k++) {
        for (l=0;l<m;l++) {
          for (q=0;q<o;q++) {

            /* Calculating zeros in y(x) and y(Q^2) */
            yx=cos(pi*(k+0.5)/n); 
            yQ=cos(pi*(l+0.5)/m);
            yr=cos(pi*(q+0.5)/o);

            /* Calculating zeros in x and Q^2 */
            zerox = pow(xmi,(0.5*(1-yx)));
            zeroQ = Qma*pow((Qma/Qmi),(0.5*(yQ-1)));
            zeror = rma*pow((rma/rmi),(0.5*(yr-1)));
            //zeror = 0.5*(rma-rmi)*yr+0.5*(rma+rmi);

            /* Function values at zeros */
            f[k][l][q]=func(zerox,zeroQ,zeror);
	    //printf("%e %e %e %e\n",zerox,zeroQ,zeror,func(zerox,zeroQ,zeror));
	    //printf("%d %d %e %e %e\n",k, l, zerox,zeroQ, f[k][l][q]);
	  }
        }
    }


    /* Calculation of Chebeshev coefficients  */
    for (i=0;i<n;i++) {
      for (j=0;j<m;j++) {
        for (q=0;q<o;q++) {
            sum = 0.0;
            for(k=0;k<n;k++) {
              for(l=0;l<m;l++)
                for(r=0;r<o;r++) {
                    sum +=f[k][l][r]*cos(pi*i*(k+0.5)/n)*cos(pi*j*(l+0.5)/m)*
		                     cos(pi*q*(r+0.5)/o);
              }
            }
        *(c+i*(m*o)+j*(o)+q) = fac*sum; 
        //printf("%d %d %d %d %e\n",i,j,q, i*(m*o)+j*(o)+q, sum);
	}
      }
    }

    //for (i=0;i<n*m*o;i++) {
    //  printf("%e\n",*(c+i));
    //}

}

/*******************************************************************************
* Chebyshev evaluation - returns the value of approximated function
*******************************************************************************/
double chebev3 (double xmi, double xma, 
               double Qmi, double Qma, 
               double rmi, double rma, 
               int n  , int m, int o,
               double* c , double x, double Q, double r) {

    double yxcheb[n];
    double yQcheb[m];
    double yrcheb[o];


    /* Transformation coefficients */
    double Ax = -2/log(xmi);
    double Bx =  1;
    double AQ =  2/log(Qma/Qmi);
    double BQ =  1-2/log(Qma/Qmi)*log(Qma);
    double Ar =  2/log(rma/rmi);
    double Br =  1-2/log(rma/rmi)*log(rma);
    //double Ar =  2/(rma-rmi);
    //double Br =  -(rma+rmi)/(rma-rmi);
     
    double yx,yQ,yr;
    double sumx = 0.0, sumQ = 0.0, sumr = 0.0;
    double sumxQ = 0.0, sumxr = 0.0, sumQr = 0.0;
    double sumxQr = 0.0;
    int i,j,q,k; /* Counters */

    /* Logarithmic transformation*/
    yx = (Ax*log(x)+Bx);
    yQ = (AQ*log(Q)+BQ);
    yr = (Ar*log(r)+Br);
    //yr = Ar*r+Br;

    yxcheb[0] = 1;
    yxcheb[1] = yx;
    yQcheb[0] = 1;
    yQcheb[1] = yQ;
    yrcheb[0] = 1;
    yrcheb[1] = yr;

    
    for (k=2;k<n;k++)
    	yxcheb[k]=2*yx*yxcheb[k-1]-yxcheb[k-2];
    for (k=2;k<m;k++)
    	yQcheb[k]=2*yQ*yQcheb[k-1]-yQcheb[k-2];
    for (k=2;k<o;k++)
    	yrcheb[k]=2*yr*yrcheb[k-1]-yrcheb[k-2];

    //*(c+i*(m*o)+j*(o)+q)

    for (i=1;i<n;i++)
        sumx += *(c+i*(m*o))*yxcheb[i];

    for (j=1;j<m;j++)
        sumQ += *(c+j*o)*yQcheb[j];

    for (q=1;q<o;q++)
        sumr += *(c+q)*yrcheb[q];

    for (i=1;i<n;i++) 
        for (j=1;j<m;j++) 
            sumxQ+=*(c+i*(m*o)+j*(o))*yxcheb[i]*yQcheb[j];

    for (i=1;i<n;i++) 
        for (q=1;q<o;q++) 
            sumxr+=*(c+i*(m*o)+q)*yxcheb[i]*yrcheb[q];

    for (j=1;j<m;j++) 
        for (q=1;q<o;q++) 
            sumQr+=*(c+j*(o)+q)*yQcheb[j]*yrcheb[q];

    for (i=1;i<n;i++) 
        for (j=1;j<m;j++) 
            for (q=1;q<o;q++) 
               sumxQr+=*(c+i*(m*o)+j*(o)+q)*yxcheb[i]*yQcheb[j]*yrcheb[q];
    //printf("# %e %e %e %e %e %e %e \n",sumx, sumQ, sumr, sumxQ, sumxr, sumQr,
    //sumxQr);
        
    return c[0]/8 + sumx/4 + sumQ/4 + sumr/4 +
           sumxQ/2 +sumxr/2 +sumQr/2 + sumxQr;
}
