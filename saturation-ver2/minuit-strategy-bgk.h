int RUN_MINUIT(void(*fcn)(int* , double*, double*, double *,unsigned*,void (*)(void) ) ){
	char command[100];
	int error_flag, istat, nvpar, npar;
	double val, edm, up;
	
	MNCOMD(*fcn,"SET STRATEGY 0",error_flag,0);
	MNCOMD(*fcn ,  "SET PRINTOUT 3",error_flag,0);	

	N_SIMPS=(int)(N_SIMPS_R*3.0/5.0);
	N_CHEB=(int)(N_CHEB_R*3.0/5.0);
	SIGMA_PREC=DGAUSS_PREC*10;
	
	MNCOMD(*fcn, "SIMPLEX",error_flag,0);
	
	
	MNCOMD(*fcn,"SET LIMITS",error_flag,0);

	MNCOMD(*fcn,"MIGRAD",error_flag,0);

	N_SIMPS=N_SIMPS_R;
	N_CHEB=N_CHEB_R;
	SIGMA_PREC=DGAUSS_PREC;

	MNCOMD(*fcn,"SET STRATEGY 1",error_flag,0);

	MNCOMD(*fcn,"MIGRAD",error_flag,0);
	
	MNSTAT(val,edm, up, nvpar,npar,istat);
	return(istat);
}
