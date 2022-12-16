double psisq_f (double *Z,void* par)  {
	double* param=(double*) par;
	
	double z=*Z;
	//double charge_sum;
	double r=param[0];
	double Q2=param[1];
	double mass2=param[2];

	double	value;
	double	z_bar =  z*z+(1-z)*(1-z);
	//double     y_bar =  (y*y)/(1+(1-y)*(1-y));
	double	Qsq_bar =  z*(1-z)*Q2+mass2;
	double	Qsq2 =  sqrt(Qsq_bar)*r;
	
	//pow(r,2) is to suppress singularity at r=0, it is compensated by the sigma
	//printf("ep=%.5e\n",Qsq2);
	
	if(Qsq2<1.0e-5){//small er approximation
		//printf("small ep\n");
		value =   (z_bar + ( mass2+ pow(2*z*(1-z),2)* Q2 )*pow(r* log(Qsq2),2) );
		
	}else{
		double	bessel_k0_2 = pow(std::cyl_bessel_k(0,Qsq2),2);
		double	bessel_k1_2 = pow(std::cyl_bessel_k(1,Qsq2),2);
		value = pow(r,2) * (z_bar * Qsq_bar * bessel_k1_2 + ( mass2 + pow(2*z*(1-z),2)* Q2 ) * bessel_k0_2);
	}
	return ( (3*value)/(2*PI*PI) );
}


Clenshaw_Curtis photon_z_integrator(32);
//Chebyshev_Gauss photon_z_integrator2(32);
double psisq_z_int(double r,double q2,double mf2){
	double res=0;
	double zmin=1.0e-15;
	double zmax = 0.5;//because psi is symmetric in z<->1-z!!
	double param[]={ r,q2,mf2};
	photon_z_integrator.name="Psi";
	//res=photon_z_integrator2.integrate(&psisq_f,(void*)param,0,zmin,1,INT_PREC/10,INT_PREC/100);
	if(zmin<1.0e-5){
		res+=dgauss(&psisq_f,(void*)param,zmin,1.0e-5,INT_PREC/10,INT_PREC/100);
	}
	if(zmax>1.0e-5){
		res+=photon_z_integrator.integrate(&psisq_f,(void*)param,1.0e-5,zmax,5,INT_PREC/10,INT_PREC/100);
	}
	
	//res=dclenshaw(&psisq_f,(void*)param,zmin,zmax, 1.0e-5/*DGAUSS_PREC*/);
	return(2.0* res);

}



