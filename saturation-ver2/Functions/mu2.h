/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////Something for the parametrizations of mu2/////////////////////////////////
///////////////////////////////Jacobian dmu2/dr and ddmu2/drdr/////////////////////////////////////////

#if (STAR==0)
/////////////////////////////////////Type 0//////////////////////////////////////////////
int compute_mu2(double r, const double * sudpar, double * mu2_arr,int opt){
	double C=sudpar[0];
	double rmax=sudpar[1];
	double r2=r*r;
	double rm2=rmax*rmax;
	
	*(mu2_arr)=C*(1.0/r2+1.0/rm2);
	
	if(opt==1){
		if( (mu2_arr[0]<LQCD2)){
			return 1;		
		}
		if( isnan(mu2_arr[0])!=0){
			return 2;		
		}
				
		return 0;
	}
	
	//first and second derivatives wrt r 
	mu2_arr[1] = -2.0*C/(r*r2);
	mu2_arr[2] = 6.0*C/(r2*r2);
	
	
	for(int i=1;i<3;i++){
		if( (isnan(mu2_arr[i])!=0) ){
			return 2;
		}
		if( (mu2_arr[i]*(pow(-1,i))<0) ){
			return 1;
		}
	}
	
	if(opt==3){
		return 0;
	}else{
		printf("opt is 1 or 3 : %d\n", opt); 
		return 3;
	}	
}


/////////////////////////////////////////////
int rmin2(double Q2 , const double* sudpar,double *rmin_2){
	double C=sudpar[0];
	double rmax=sudpar[1];
	double frac=Q2*rmax*rmax/C;
	if(C<0){
		printf("negative C, %.3e",C);
		return 1;
	}
	if(frac<1){
		return 9;
	}
	
	
	*rmin_2=rmax*rmax/(frac-1);
	
	if(rmin2<0){
		printf("r min negative");
		getchar();
		return 1;
	}
	return 0;
}


#elif (STAR==1)

///////////////////////////////////////type1////////////////////////////////////////////

int compute_mu2(double r, const double * sudpar, double * mu2_arr, int opt){
	double C=sudpar[0];
	double rmax=sudpar[1];
	double r2=r*r;
	double rm2=rmax*rmax;
	
	double exprrmax=exp(-r2/rm2);
	
	mu2_arr[0]=C/(rm2*(1.0-exprrmax ));
	
	if(opt==1){
		if( (mu2_arr[0]<LQCD2)){
			return 1;		
		}
		if( isnan(mu2_arr[0])!=0){
			return 2;		
		}
		return 0;
	}
	
	double jac=-2.0*C/(rm2*rm2)*( exprrmax/((1.0-exprrmax)*(1.0-exprrmax))); 
	//first and second derivatives wrt r 
	mu2_arr[1] = jac*r;
	mu2_arr[2] = jac*( 1.0- 2* (r2/rm2) *((1.0+exprrmax)/(1.0-exprrmax)) ) ;	
	
	for(int i=1;i<3;i++){
		if( (isnan(mu2_arr[i])!=0) ){
			return 2;
		}
		if( (mu2_arr[i]*(pow(-1,i))<0) ){
			return 1;
		}
	}
	if(opt==3){
		return 0;
	}else{
		printf("opt is 1 or 3 : %d\n", opt); 
		return 3;
	}
}


/////////////////////////////////////////////
int rmin2(double Q2 ,const double* sudpar, double *rmin_2){
	double C=sudpar[0];
	double rmax=sudpar[1];
	double frac=C/(Q2*rmax*rmax); 
	if(C<0){
		printf("negative C, %.3e",C);
		return 1;
	}
	if((frac>1)){
		return 9;	
	}
	
	*rmin_2=-rmax*rmax*log( 1.0-frac);
	
	if(rmin2<0){
		printf("r min negative");
		getchar();
		return 1;
	}
	return 0;
}
#endif

