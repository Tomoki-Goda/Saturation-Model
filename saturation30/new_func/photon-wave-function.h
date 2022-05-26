//#include<stdio.h>
//#include<math.h>
//#include"constants.h"
//


extern double dbesk0_(double*);
extern double dbesk1_(double*);

double psisq_f (double r, double z, double Q2, unsigned char flavourtype/*, unsigned dataform */)  {
//	double z=exp(-z_mod/(1-z_mod));
//	double jac=z_mod/pow(1-z,2);		
	double charge_sum;
	double mass2;

	switch(flavourtype){
		case 'l':
			charge_sum=5.0/6.0;
			mass2=MASS_L2;
			break;
		case 's':
			charge_sum=1.0/6.0;
			mass2=MASS_S2;
			break;
		case 'c':
			charge_sum=2.0/3.0;
			mass2=MASS_C2;
			break;
		case 'b':
			charge_sum=1.0/6.0;
			mass2=MASS_B2;
			break;
		default:
			printf("wrong input %c\n",flavourtype);
			charge_sum=5.0/6.0;
			mass2=MASS_L2;
	}

	double	z_bar =  z*z+(1-z)*(1-z);
	//double     y_bar =  (y*y)/(1+(1-y)*(1-y));
	double	Qsq_bar =  z*(1-z)*Q2+mass2;
	double	Qsq2 =  sqrt(Qsq_bar)*r;
//	double	norm_f;
	//external bessels functions, declared  elsewhere
	double	bessel_k0 = dbesk0_(&Qsq2);
	double	bessel_k1 = dbesk1_(&Qsq2);
	double	value;
			
// F_2 form
	//pow(r,2) is to suppress singularity at r=0, it is compensated by the sigma
	//Q2 comes from the factor to multiply F2 to get cross-section.
	value = pow(r,3)*Q2*  (charge_sum)*NORM*(z_bar*Qsq_bar*bessel_k1*bessel_k1 + ( mass2 +4*Q2*z*z*(1-z)*(1-z))*bessel_k0*bessel_k0);
	return (value);
}
double psisq_z_integrand(double z,double **par){
	const double r=(*(*par));
	const double Q2=(*((*par)+1));
	unsigned char f ;//(*(par+1))
	switch((int)(*((*par)+2))){
		//just anything to convert double to char 
		case 0:
			f='l';
			break;
		case 1:
			f='s';
			break;
		case 2:
			f='c';
			break;
		case 3:
			f='b';
			break;
		default:
		printf("wrong input %f\n",*((*par)+2));
		f='l';
	}
	
	double val=psisq_f(r,z,Q2,f);
	//printf( " r %f  z %f  q2 %f ->%f \n", r,z,Q2,val);
	return(val);
}

double psisq_z_int(double r,double Q2,unsigned char f){
	double f_n;
	switch(f){//.2 for safety when converting from double to char
                case 'l':
                        f_n=0.2;
                        break;
                case 's':
                        f_n=1.2;
                        break;
                case 'c':
                        f_n=2.2;
                        break;
                case 'b':
                        f_n=3.2;
                        break;
                default:
			printf("wrong input %c\n",f);
			f_n=0.2;
        }	
	double param[3];//={r,Q2,f_n};
	*(param)=r;
	*(param+1)=Q2;
	*(param+2)=f_n;
	
	double *par[1];// ={param};
	*(par)=param;
	double res=0;
	double ep =1.0e-6;
	//printf("r:%f %f, Q: %f %f, %f %c\n",(*(*par)),r,(*((*par)+1)) ,Q2,(*((*par)+2)),f);
	simpson1d(&psisq_z_integrand,par,0.0+ep,1.0-ep,&res);//1.0e-15,1.0-1.0e-15/*1.0-1.0e-10*/,&res);
	//printf("Q2: %f, r: %f, res: %f  %c\n",Q2,r,res,f);
	return(res);

}

