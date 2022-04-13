//
//     Author: Sebastian Sapeta
//       Date: August 2005
//

#define  MX 22   /* The order of polynomial in x variable */
#define  MQ 20   /* The order of polynomial in Q^2 variable */

double coef[MX][MQ];     /* The coefficients vector */
double xmin =  1.0e-07;  /* x lower limit */
double xmax =  1.0e+00;  /* x upper limit */
double Qmin =  1.0e-01;  /* Q^2 upper limit */
double Qmax =  1.0e+03;  /* Q^2 upper limit */

//#define  NX 20   /* The order of polynomial in x variable */
//#define  NR 20   /* The order of polynomial in x variable */
//double rmin =  1.0e-06;  /* r lower limit */
//double rmax =  3.0e+01;  /* r upper limit */
//double coef2[NX][NR];     /* The coefficients vector */

/*******************************************************************************
*  MAX function - returns max value of two integers
*******************************************************************************/
int max (int a, int b) {

    int i;
    if(a>b) 
        i = 1;
    else if (b>a)
        i = 0;
    return i*a + (1-i)*b;
}


/*******************************************************************************
* Chebyshev - coefficients calculation
*******************************************************************************/
extern void chebft (double xmi, double xma, 
             double Qmi, double Qma, 
             int n, int m,
             double* c , double (*func)(double,double)) {
    int    k,l,i,j;  /* Counters */
    double fac = 4.0/(n*m);
    double sum;
    double f[n][m]; /* Function values array */ 

    double yx; /* Polynomial's zeros in y(x)*/
    double yQ; /* Polynomial's zeros in y(Q^2) */

    double zerox; /* Polynomial's zeros in x */
    double zeroQ; /* Polynomial's zeros in Q^2 */



    /* Calculation of function values at zeros */
    for (k=0;k<n;k++) {
        for (l=0;l<m;l++) {

            /* Calculating zeros in y(x) and y(Q^2) */
            yx=cos(pi*(k+0.5)/n); 
            yQ=cos(pi*(l+0.5)/m);

            /* Calculating zeros in x and Q^2 */
            zerox = pow(xmi,(0.5*(1-yx)));
            zeroQ = Qma*pow((Qma/Qmi),(0.5*(yQ-1)));
            //zeroQ = 0.5*(Qma-Qmi)*yQ+0.5*(Qma+Qmi);

            /* Function values at zeros */
            f[k][l]=func(zerox,zeroQ);  
	    //printf("%d %d %e %e %e\n",k, l, zerox,zeroQ, f[k][l]);
        }
    }


    /* Calculation of Chebeshev coefficients  */
    for (i=0;i<n;i++) {
        for (j=0;j<m;j++) {
            sum = 0.0;
            for(k=0;k<n;k++) {
                for(l=0;l<m;l++)
                    sum +=f[k][l]*cos(pi*i*(k+0.5)/n)*cos(pi*j*(l+0.5)/m);
            }
        c[i+j*m]=fac*sum; 
        }
    }

}

/*******************************************************************************
* Chebyshev Polynomials - returns the value of i-th polynomial
*******************************************************************************/
double chebpol (int i, double x) {

    double T[max(MX,MQ)]; /* Array of polynomials' values */
    int k; /* Counter */
  
    T[0]=1;
    T[1]=x;
  
    for(k=2;k<=i;k++)
        T[k]=2*x*T[k-1]-T[k-2];

    return T[i];
}

/*******************************************************************************
* Chebyshev evaluation - returns the value of approximated function
*******************************************************************************/
double chebev (double xmi, double xma, 
               double Qmi, double Qma, 
               int n  , int m,
               double* c , double x, double Q) {

    double yxcheb[n];
    double yQcheb[m];


    /* Transformation coefficients */
    double Ax = -2/log(xmi);
    double Bx =  1;
    double AQ =  2/log(Qma/Qmi);
    double BQ =  1-2/log(Qma/Qmi)*log(Qma);
    //double AQ = 2.0/(Qma-Qmi);
    //double BQ = -(Qma+Qmi)/(Qma-Qmi);
     
    double yx,yQ;
    double sumx = 0.0, sumQ = 0.0, sumxQ = 0.0; /* Partial summs */
    int i,j,k; /* Counters */

    /* Logarithmic transformation*/
    yx = (Ax*log(x)+Bx);
    yQ = (AQ*log(Q)+BQ);
    //yQ = AQ*Q+BQ;

    yxcheb[0] = 1;
    yxcheb[1] = yx;
    yQcheb[0] = 1;
    yQcheb[1] = yQ;

    
    for (k=2;k<n;k++) 
    	yxcheb[k]=2*yx*yxcheb[k-1]-yxcheb[k-2];   
    for (k=2;k<m;k++) 
    	yQcheb[k]=2*yQ*yQcheb[k-1]-yQcheb[k-2];   

    for (i=1;i<n;i++) 
        sumx += c[i]*yxcheb[i];

    for (j=1;j<m;j++) 
        sumQ += c[j*m]*yQcheb[j];

    for (i=1;i<n;i++) 
        for (j=1;j<m;j++) 
            sumxQ+=c[i+j*m]*yxcheb[i]*yQcheb[j];
        
    return 0.25*c[0]+0.5*sumx+0.5*sumQ+sumxQ;
}
/*******************************************************************************
* Chebyshev evaluation - returns the value of approximated function
*******************************************************************************/
double chebev2 (double xmi, double xma, 
               double Qmi, double Qma, 
               int n      , int m     ,
               double c[MX][MQ] , double x, double Q) {

    /* Transformation coefficients */
    double Ax = -2/log(xmi);
    double Bx =  1;
    double AQ =  2/log(Qma/Qmi);
    double BQ =  1-2/log(Qma/Qmi)*log(Qma);
     
    double yx,yQ;
    double sumx = 0.0, sumQ = 0.0, sumxQ = 0.0; /* Partial summs */
    int i,j; /* Counters */

    /* Logarithmic transformation*/
    yx = (Ax*log(x)+Bx);
    yQ = (AQ*log(Q)+BQ);

    for (i=1;i<n;i++) 
        sumx+=c[i][0]*chebpol(i,yx);

    for (j=1;j<m;j++) 
        sumQ += c[0][j]*chebpol(j,yQ);

    for (i=1;i<n;i++) 
        for (j=1;j<m;j++) 
            sumxQ+=c[i][j]*chebpol(i,yx)*chebpol(j,yQ);
        
    return 0.25*c[0][0]+0.5*sumx+0.5*sumQ+sumxQ;
}


