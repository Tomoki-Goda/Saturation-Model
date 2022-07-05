/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////Something for the parametrizations of mu2/////////////////////////////////
///////////////////////////////Jacobian dmu2/dr and ddmu2/drdr/////////////////////////////////////////

#if STAR ==0
/////////////////////////////////////Type 0//////////////////////////////////////////////
int compute_mu2(double r, const double * sudpar, double * const mu2_arr,int opt){
	double C=sudpar[0];
	double rmax=sudpar[1];
	double r2=r*r;
	double rm2=rmax*rmax;
	mu2_arr[0]=C*(1.0/r2+1.0/rm2);
	if(opt==1){
		return 0;
	}
	mu2_arr[1] = -2.0*C/(r*r2);
	mu2_arr[2] = 6.0*C/(r2*r2);
	return 0;	
}


/////////////////////////////////////////////
double rmin2(double Q2 ,const double* sudpar){
	double C=sudpar[0];
	double rmax=sudpar[1];
	
	double rmin2=-rmax*rmax*C/(C-Q2*rmax*rmax);
	return rmin2;
}


#elif STAR ==1

///////////////////////////////////////type1////////////////////////////////////////////

int compute_mu2(double r, const double * sudpar, double * const mu2_arr, int opt){
	double C=sudpar[0];
	double rmax=sudpar[1];
	double r2=r*r;
	double rm2=rmax*rmax;
	
	double exprrmax=exp(-r2/rm2);
	
	mu2_arr[0]=C/(rm2*(1.0-exprrmax ));
	if(opt==1){
		return 0;
	}
	double jac=-2.0*C/(rm2*rm2)*( exprrmax/pow(1.0-exprrmax,2)); 
	
	mu2_arr[1] = jac*r;
	mu2_arr[2] = jac*( 1.0- 2* (r2/rm2) *((1.0+exprrmax)/(1.0-exprrmax)) ) ;	
	return 0;
}


/////////////////////////////////////////////
double rmin2(double Q2 ,const double* sudpar){
	double C=sudpar[0];
	double rmax=sudpar[1];
	
	double rmin2=-rmax*rmax*log( 1.0-C/(Q2*rmax*rmax) );
	return rmin2;
}
#endif

