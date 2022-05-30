//#include<stdio.h>
//#include<math.h>
//#include"constants.h"
//


extern double dbesk0_(double*);
extern double dbesk1_(double*);

extern double dgquad_(double (*)(double*), double*,double*,int*  );

static double Q2;
static double R;
static unsigned flavourtype;

double psisq_f (double *Z)  {
//double psisq_f (double r, double z, double Q2, unsigned char flavourtype/*, unsigned dataform */)  {
//	double z=exp(-z_mod/(1-z_mod));
//	double jac=z_mod/pow(1-z,2);	
	double z=*Z;	
	double charge_sum;
	double mass2;

	switch(flavourtype){
		case 0:
			charge_sum=5.0/6.0;
			mass2=MASS_L2;
			break;
		case 1:
			charge_sum=1.0/6.0;
			mass2=MASS_S2;
			break;
		case 2:
			charge_sum=2.0/3.0;
			mass2=MASS_C2;
			break;
		case 3:
			charge_sum=1.0/6.0;
			mass2=MASS_B2;
			break;
		default:
			printf("psisq_f::wrong input %d\n",flavourtype);
			charge_sum=5.0/6.0;
			mass2=MASS_L2;
	}

	double	z_bar =  z*z+(1-z)*(1-z);
	//double     y_bar =  (y*y)/(1+(1-y)*(1-y));
	double	Qsq_bar =  z*(1-z)*Q2+mass2;
	double	Qsq2 =  sqrt(Qsq_bar)*R;
//	double	norm_f;
	//external bessels functions, declared  elsewhere
	double	bessel_k0_2 = pow(dbesk0_(&Qsq2),2);
	double	bessel_k1_2 = pow(dbesk1_(&Qsq2),2);
	//double	bessel_k0_reg = 1/r    - dbesk0_(&Qsq2);
	//double	bessel_k1_reg = 1/Qsq2 - dbesk1_(&Qsq2);
	
	double	value;
			
// F_2 form
	//pow(r,2) is to suppress singularity at r=0, it is compensated by the sigma
	//Q2 comes from the factor to multiply F2 to get cross-section.
	value = pow(R,2) * Q2 *  (charge_sum) * NORM * (z_bar * Qsq_bar * bessel_k1_2 + ( mass2 + pow(2*z*(1-z),2)* Q2 ) * bessel_k0_2);
	//value = Q2 *  (charge_sum) * NORM * (z_bar * pow(1-Qsq2*bessel_k1_reg,2) + ( mass2 + pow(2*z*(1-z),2)* Q2 ) * pow(1-R*bessel_k0_reg,2));
	return (value);
}



double psisq_z_int(double r,double q2,unsigned f){
	flavourtype=f;
	R=r;
	Q2=q2;
	double res=0;
	double zmin=1.0e-5;
	double zmax =1.0 - 1.0e-5;
	int N=96;//64;
	
	res=dgquad_(&psisq_f,&zmin,&zmax,&N);
	//printf("%f\n",res);
	return(res);

}




