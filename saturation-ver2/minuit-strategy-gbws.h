int RUN_MINUIT(void(*fcn)(int* , double*, double*, double *,unsigned*,void (*)(void) ) ){
	char command[100];
	int error_flag, istat, nvpar, npar;
	double val, edm, up;
	
	MNCOMD(*fcn,"SET STRATEGY 0",error_flag,0);
	MNCOMD(*fcn ,  "SET PRINTOUT 2",error_flag,0);

	N_SIMPS=(int)((N_SIMPS_R)/4);
	N_CHEB=(int)((N_CHEB_R)/4);
	SIGMA_PREC=DGAUSS_PREC*50;

	MNCOMD(*fcn, "SIMPLEX",error_flag,0);

	N_SIMPS=(int)((N_SIMPS_R)/2);
	N_CHEB=(int)((N_CHEB_R)/2);
	SIGMA_PREC=DGAUSS_PREC*10;

	MNCOMD(*fcn, "SIMPLEX 150 1.0D0",error_flag,0);
	
	
	MNCOMD(*fcn,"FIX 5",error_flag,0);
	
	MNCOMD(*fcn,"SET LIMITS",error_flag,0);
	
	

	MNCOMD(*fcn,"MIGRAD",error_flag,0);

	N_SIMPS=N_SIMPS_R;
	N_CHEB=N_CHEB_R;
	SIGMA_PREC=DGAUSS_PREC;
	
	MNCOMD(*fcn,"SET STRATEGY 1",error_flag,0);
	

	MNCOMD(*fcn,"FIX 4",error_flag,0);
	MNCOMD(*fcn,"FIX 5",error_flag,0);
	
	MNCOMD(*fcn,"HESSE ",error_flag,0);
	MNCOMD(*fcn,"MIGRAD",error_flag,0);
	
	MNCOMD(*fcn,"RELEASE 4",error_flag,0);
	
	MNCOMD(*fcn,"HESSE ",error_flag,0);
	MNCOMD(*fcn,"MIGRAD",error_flag,0);
	
	MNCOMD(*fcn,"RELEASE 5",error_flag,0);
	MNCOMD(*fcn,"SET LIMITS 5 1.0D0 20.0D0",error_flag,0);

	for(int i=0 ; i<3;i++){
		MNCOMD(*fcn,"HESSE ",error_flag,0);
		
		MNCOMD(*fcn,"MIGRAD",error_flag,0);
		MNSTAT(val,edm, up, nvpar,npar,istat);
		if(istat==3){
			break;
		}else{	
			MNCOMD(*fcn,"SIMPLEX 150 0.1D-3",error_flag,0);
			if(i==2){
				MNCOMD(*fcn,"HESSE ",error_flag,0);
				MNCOMD(*fcn,"MIGRAD",error_flag,0);
				MNSTAT(val,edm, up, nvpar,npar,istat);
			}
		}
	}

	
	MNSTAT(val,edm, up, nvpar,npar,istat);
	return(istat);
}
