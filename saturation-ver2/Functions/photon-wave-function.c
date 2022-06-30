#include<stdio.h>
#include<math.h>


//#include"./control_tmp.h"
#include"control.h"
#include"control-default.h"
#include"constants.h"

//#include"simpson-integral.h"


extern double dbesk0_(double*);
extern double dbesk1_(double*);

extern double dgquad_(double (*)(double*), double*,double*,int*  );

extern void simpson1dA(double(*)(double ,double**),double**,double,double,int,double*,double*);


#if SIMPS_Z_INT==0
static double Q2;
static double R;
static unsigned flavourtype;

double psisq_f (double *Z)  {
	double z=*Z;
#else
double psisq_f (double R, double z, double Q2, unsigned flavourtype/*, unsigned dataform */)  {
	
#endif
//	double z=exp(-z_mod/(1-z_mod));
//	double jac=z_mod/pow(1-z,2);	
		
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
			getchar();
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


#if SIMPS_Z_INT==0
double psisq_z_int(double r,double q2,unsigned f){
	//printf("%d\n",f);
	flavourtype=f;
	R=r;
	Q2=q2;
	double res=0;
	double zmin=1.0e-5;
	//double zmax =1.0 - 1.0e-5;
	double zmax = 0.5;//because psi is symmetric in z<->1-z!!
	int N=96;//64;
	
	res=dgquad_(&psisq_f,&zmin,&zmax,&N);
	//printf("%f\n",res);
	return(2.0* res);

}
#else
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
double psisq_z_integrand(double z,double **par){
	const double r=(*(*par));
	const double Q2=(*((*par)+1));
	unsigned f =(unsigned) (*((*par)+2) +0.5);//(*(par+1))
	
	double val=psisq_f(r,z,Q2,f);
	//printf( " r %f  z %f  q2 %f ->%f \n", r,z,Q2,val);
	return(val);
}

double psisq_z_int(double r,double Q2,unsigned f){
		
	double param[3];
	*(param)=r;
	*(param+1)=Q2;
	*(param+2)=(double)f;
	
	double *par[1];
	*(par)=param;
	double res=0;
	double ep =1.0e-6;
	double err=0;
	simpson1dA(&psisq_z_integrand,par,0.0+ep,0.5,100,&res,&err);//zmax is 0.5 because psi is symmetric in z<->1-z!!
	
	return(2.0* res);

}


#endif

