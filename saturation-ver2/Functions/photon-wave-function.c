#include<stdio.h>
#include<math.h>


//#include"./control_tmp.h"
#include"control.h"
#include"control-default.h"
#include"constants.h"

//#include"simpson-integral.h"


extern double dbesk0_(double*);
extern double dbesk1_(double*);

extern double dgquad_(double (*)(const double*), double*,double*,int*  );
extern double dgauss_(double (*)(const double*), double*,double*,double *  );
extern void simpson1dA(double(*)(double ,double**),double**,double,double,int,double*,double*);
extern double dadapt_(double(* )(const double*),double*,double*,int*,double*,double* ,double*,double*);

#if SIMPS_Z_INT==0
static double Q2;
static double R;
static unsigned flavourtype;

double psisq_f (const double *Z)  {
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
	double	value;
	double	z_bar =  z*z+(1-z)*(1-z);
	//double     y_bar =  (y*y)/(1+(1-y)*(1-y));
	double	Qsq_bar =  z*(1-z)*Q2+mass2;
	double	Qsq2 =  sqrt(Qsq_bar)*R;
	
	//pow(r,2) is to suppress singularity at r=0, it is compensated by the sigma
	//printf("ep=%.5e\n",Qsq2);
	if(Qsq2<1.0e-5){//small er approximation
		//printf("small ep\n");
		value =   (z_bar + ( mass2+ pow(2*z*(1-z),2)* Q2 )*pow(R* log(Qsq2),2) );
		
	}else{
		double	bessel_k0_2 = pow(dbesk0_(&Qsq2),2);
		double	bessel_k1_2 = pow(dbesk1_(&Qsq2),2);
		value = pow(R,2) * (z_bar * Qsq_bar * bessel_k1_2 + ( mass2 + pow(2*z*(1-z),2)* Q2 ) * bessel_k0_2);
	}
	//Q2 comes from the factor to multiply F2 to get cross-section.
	return (Q2 *  (charge_sum) * NORM *value);
}


#if SIMPS_Z_INT==0
double psisq_z_int(double r,double q2,unsigned f){
	//printf("%d\n",f);
	flavourtype=f;
	R=r;
	Q2=q2;
	double res=0;
	double zmin=0.0;
	//double zmin=1.0e-10;
	//double zmax =1.0 - 1.0e-5;
	double zmax = 0.5;//because psi is symmetric in z<->1-z!!

	
	//int N=96;
	//double step=(zmax-zmin)/5;
	
	//res=0;
	//for(int i=0;i<5;i++){
	//	zmax=zmin+step;
	//	res+=dgquad_(&psisq_f,&zmin,&zmax,&N);
	//	zmin=zmax;
	//}
	//res=0;
	//double low, high;
	//for(int i=0;i<5;i++){
	//	low=zmin*pow((zmax/zmin),((double)(i))/5);
	//	high=zmin*pow((zmax/zmin),((double)(i+1))/5);
	//	res+=dgquad_(&psisq_f,&low,&high,&N);
	//}
	//res=dgquad_(&psisq_f,&zmin,&zmax,&N);
	//////////////////////////////////
	//double N=DGAUSS_PREC/100;
	//res=dgauss_(&psisq_f,&zmin,&zmax,&N);
	//double step=(zmax-zmin)/1;
	res=0;
	//double low=zmin, high;
	//for(int i=0;i<1;i++){
	//	high=low+step;
		//low=zmin*pow((zmax/zmin),((double)(i))/5);
		//high=zmin*pow((zmax/zmin),((double)(i+1))/5);
	//	res+=dgauss_(&psisq_f,&low,&high,&N);
	//	low=high;
	//}
	
	////////////////////////////////////
	int seg=5;
	double NRel=DGAUSS_PREC;
	double NAbs=0;
	double error=0;
	dadapt_(&psisq_f,&zmin,&zmax,&seg ,&NRel, &NAbs, &res, &error)	;
	/////////////////////////////////////////
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
	printf("psisq_z_int no longer supported \n");
	getchar();
		
	double param[3];
	*(param)=r;
	*(param+1)=Q2;
	*(param+2)=(double)f;
	
	double *par[1];
	*(par)=param;
	double res=0;
	double ep =1.0e-10;
	double err=0;
	simpson1dA(&psisq_z_integrand,par,0.0+ep,0.5,100,&res,&err);//zmax is 0.5 because psi is symmetric in z<->1-z!!
	
	return(2.0* res);

}


#endif

