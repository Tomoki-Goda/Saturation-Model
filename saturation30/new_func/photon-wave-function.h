//#include<stdio.h>
//#include<math.h>
//#include"constants.h"
//


extern double dbesk0_(double*);
extern double dbesk1_(double*);

double psisq_f (double r, double z, double Q2, unsigned char flavourtype/*, unsigned dataform */)  {
			
	double charge_sum;
	double mass;

	switch(flavourtype){
		case 'l':
			charge_sum=5.0/6.0;
			mass=MASS_L2;
			break;
		case 's':
			charge_sum=1.0/6.0;
			mass=MASS_S2;
			break;
		case 'c':
			charge_sum=2.0/3.0;
			mass=MASS_C2;
			break;
			case 'b':
			charge_sum=1.0/6.0;
			mass=MASS_B2;
			break;
	}

	double	z_bar =  z*z+(1-z)*(1-z);
	//double     y_bar =  (y*y)/(1+(1-y)*(1-y));
	double	Qsq_bar =  z*(1-z)*Q2+mass;
	double	Qsq2 =  sqrt(Qsq_bar)*r;
	double	norm_f;
	//external bessels functions, declared  elsewhere
	//printf("%f,  %f,  ",r,Qsq2);
	double	bessel_k0 = dbesk0_(&Qsq2);
	double	bessel_k1 = dbesk1_(&Qsq2);
	double	value;
			
	switch (1/*dataform*/) {
		case 1: // F_2 form
			value = (charge_sum)*NORM*Q2*(z_bar*Qsq_bar*bessel_k1*bessel_k1 + ( mass +4*Q2*z*z*(1-z)*(1-z))*bessel_k0*bessel_k0);
			break;
		case 2: // F_L form 
			value = (charge_sum)*NORM * Q2*4*Q2*z*z*(1-z)*(1-z)*bessel_k0*bessel_k0;
		        break;
	}
	/*
	if (photo==1) {
		norm_f = (charge_sum)*(charge_sum);
		value = norm_f*z_bar*Qsq_bar*bessel_k1*bessel_k1;
	}
	*/
	//printf("%f\n", value);
	return (value);
}

