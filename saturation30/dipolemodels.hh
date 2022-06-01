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

    double mub2=C/(pow(r,2)) + mu02;
   
    if (mu2 < Lambda2 || mub2 < Lambda2|| mu2 < mub2) {/*printf("1"); */return(0.0);}; 

    double CA = 3.0;
    double b0 = (33-2*n_f)/12;

    double val = CA/(2*b0*PI)*
                 (log(mu2/Lambda2)*log(log(mu2/Lambda2)/log(mub2/Lambda2)) - 
	         log(mu2/mub2));

    return val;
}

/******************************************************************************
 *non perturbative sudakov...
 * ***************************************************************************/
double sudakov_np(double  r,double mu2){
	double val=g1 * pow(r,2.0)/(2.0);//+ g2 * ( log(mu2/pow(Q0,2.0)) * log((pow(bmax,2.0)+pow(r,2.0))/pow(bmax,2.0))/4.0 );
	return(val);
};

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
   return exp(chebev3(xmin,xmax,Qmin,Qmax,rmin,rmax,NX,NQ,NR,*(*coef3),x,q2,r));
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
   xgpdf_cheb = chebev(xmin,xmax,Qmin,Qmax,MX,MQ,(*coef),xmod,mu2); 


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
* The dipole-proton BGK cross section modified  - Chebeshev approximation
*******************************************************************************/
double sigma_bgk_cheb_mod (double r) {

   double mu2 = C/(r*r) + mu02;
   //double mu2 = mu02/(1-exp(-mu02*r*r/C));
   double exnum,cacf; /* The numerator in the exp function */
   double xgpdf_cheb;

   xgpdf_cheb = chebev(xmin,xmax,Qmin,Qmax,MX,MQ,*coef,xmod,mu2); 
   exnum = 0.389379*pi*pi*r*r*alpha_s(mu2)*xgpdf_cheb;
  
   /* If the xgpdf is negative print the warning message */
   if (exnum <0) 
   	printf("!!! %e   %e\n", xmod, mu2);

   cacf = 9.0/4.0; 
   return sigma_0*(1-exp(-exnum*cacf/(3*sigma_0)));
}

