/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////Something for the parametrizations of mu2/////////////////////////////////
///////////////////////////////Jacobian dmu2/dr and ddmu2/drdr/////////////////////////////////////////

#if STAR ==0
/////////////////////////////////////Type 0//////////////////////////////////////////////
double rmu2( double r ,double* sudpar){
	double C=sudpar[0];
	double rmax=sudpar[1];
	
	double mu2=C*(pow(r,-2.0)+pow(rmax,-2.0));
	return mu2;
}
//////////////////////////
double rmu2_jac_first( double r ,double* sudpar){
	double C=sudpar[0];
	double rmax=sudpar[1];
	
	double	jac= -2.0*C/pow(r,3.0);
	return jac;
}
//////////////////////////
double rmu2_jac_second( double r  ,double* sudpar){
	double C=sudpar[0];
	double rmax=sudpar[1];
	
	double jac=6.0*C/pow(r,4.0);
	return jac;
}

/////////////////////////////////////////////
double rmin2(double Q2 ,double* sudpar){
	double C=sudpar[0];
	double rmax=sudpar[1];
	
	double rmin2=-rmax*rmax*C/(C-Q2*rmax*rmax);
	return rmin2;
}

#elif STAR ==1
///////////////////////////////////////type1////////////////////////////////////////////
double rmu2( double r  ,double* sudpar){
	double C=sudpar[0];
	double rmax=sudpar[1];
	
	double exprrmax=exp(-pow(r/rmax,2));
	double mu2=C/(rmax*rmax*(1.0-exprrmax ));
	return mu2;
}
//////////////////////////
double rmu2_jac_first( double r  ,double* sudpar){
	double C=sudpar[0];
	double rmax=sudpar[1];
	
	double exprrmax=exp(-pow(r/rmax,2));
	double jac=-2.0*C*(r/pow(rmax,4))*( exprrmax/pow(1.0-exprrmax,2));
	return jac;
}
//////////////////////////
double rmu2_jac_second( double r ,double* sudpar){
	double C=sudpar[0];
	double rmax=sudpar[1];
	
	double exprrmax=exp(-pow(r/rmax,2));
	double jac =- 2.0*(C/pow(rmax,4))*(exprrmax/pow(1.0-exprrmax,2))*( 1.0- 2* pow(r/rmax,2)*((1.0+exprrmax)/(1.0-exprrmax)) ) ;
	return jac;
}

/////////////////////////////////////////////
double rmin2(double Q2 ,double* sudpar){
	double C=sudpar[0];
	double rmax=sudpar[1];
	
	double rmin2=-rmax*rmax*log( 1.0-C/(Q2*rmax*rmax) );
	return rmin2;
}
#endif

