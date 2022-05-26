//#include<stdio.h>
//#include<math.h>

//#include"./dipole-cross-section.h"
//#include"./photon-wave-function.h"
//#include"./simpson-integral.h"
//#include"./constants.h"


//#define FLAVOUR 2 
#if MODEL==0
#define SIGMA sigma_gbw
#elif MODEL==1
#define SIGMA sigma_bgk
#elif MODEL==2
#define SIGMA sigma_gbs
#endif
double mod_x(double x, double Q2, unsigned char flavour) {
	double m_fsq;

	switch (flavour) {
	case 'l':
		m_fsq = MASS_L2;
		break;
	case 's':
		m_fsq = MASS_S2;
		break;
	case 'c':
		m_fsq = MASS_C2;
		break;
	case 'b':
		m_fsq = MASS_B2;
		break;
	default:
		printf("wrong input %c\n",flavour);
		m_fsq = MASS_L2;
	}
	return (x * (1.0 +( 4.0 * (m_fsq/Q2)) ));
}
/*
double sigma_integrand(double R,double z,double ** par){

	//R is r/(1+r) such that integration between (0, inf) -> (0,1).
	//double r=R/(1-R);
	double r=R;
	//double jacob=2.0*PI;
	//double jacob=1.0/pow(1-R,2);	
	double x=(*(*par));
	double Q2=(*( (*par)+1 ));
	
	double* param=( *(par+1) );
	double value=0;
	//double x_mod[4];
	
	#if FLAVOUR==0	
       		value=SIGMA(r,mod_x(x,Q2,'l'), Q2, param) * psisq_f(r, z, Q2, 'l');
                value+=SIGMA(r,mod_x(x,Q2,'s'), Q2, param) * psisq_f(r, z, Q2, 's');
	#elif FLAVOUR==1
		       	value=SIGMA(r,mod_x(x,Q2,'l'), Q2, param) * psisq_f(r, z, Q2, 'l');
                	value+=SIGMA(r,mod_x(x,Q2,'s'), Q2, param) * psisq_f(r, z, Q2, 's');
                	value+=SIGMA(r,mod_x(x,Q2,'c'), Q2, param) * psisq_f(r, z, Q2, 'c');
	#elif FLAVOUR==2
        		value=SIGMA(r,mod_x(x,Q2,'l'),Q2, param) * psisq_f(r, z, Q2, 'l');
			value+=SIGMA(r,mod_x(x,Q2,'s'), Q2, param) * psisq_f(r, z, Q2, 's');
			value+=SIGMA(r,mod_x(x,Q2,'c'), Q2, param) * psisq_f(r, z, Q2, 'c');
			value+=SIGMA(r,mod_x(x,Q2,'b'), Q2, param) * psisq_f(r, z, Q2, 'b');

//	#elif FLAVOUR==100
//			value= psisq_f(r, z, Q2, 'l');
//                     value+= psisq_f(r, z, Q2, 's');
//                     value+= psisq_f(r, z, Q2, 'c');
//                     value+= psisq_f(r, z, Q2, 'b');
	#endif

	//return(value*r*pow(1-R,-2));
	return(value*r);
}
*/

double sigma_r_integrand(double R,double ** par){

        //R is r/(1+r) such that integration between (0, inf) -> (0,1).
	double r=R/(1-R);
//	double r=R;

        double x=(*(*par));
        double Q2=(*( (*par)+1 ));

        double* param=( *(par+1) );
        double value=0;
        //double x_mod[4];
	//printf("%f\n",Q2);
        #if FLAVOUR==0  
                value=SIGMA(r,mod_x(x,Q2,'l'), Q2, param) * psisq_z_int(r, Q2, 'l');
                value+=SIGMA(r,mod_x(x,Q2,'s'), Q2, param) * psisq_z_int(r, Q2, 's');
                value*=1.0/pow(r,2); //1 from polar coordinates, -2 from compensating psi
        #elif FLAVOUR==1
                        value=SIGMA(r,mod_x(x,Q2,'l'), Q2, param) * psisq_z_int(r, Q2, 'l');
                        value+=SIGMA(r,mod_x(x,Q2,'s'), Q2, param) * psisq_z_int(r, Q2, 's');
                        value+=SIGMA(r,mod_x(x,Q2,'c'), Q2, param) * psisq_z_int(r, Q2, 'c');
                        value*=1.0/pow(r,2);
        #elif FLAVOUR==2
                        value=SIGMA(r,mod_x(x,Q2,'l'),Q2, param) * psisq_z_int(r, Q2, 'l');
                        value+=SIGMA(r,mod_x(x,Q2,'s'), Q2, param) * psisq_z_int(r, Q2, 's');
                        value+=SIGMA(r,mod_x(x,Q2,'c'), Q2, param) * psisq_z_int(r, Q2, 'c');
                        value+=SIGMA(r,mod_x(x,Q2,'b'), Q2, param) * psisq_z_int(r, Q2, 'b');
                        value*=1.0/pow(r,2);
        #endif
	
	//printf("Q2; %f, x; %f, r; %f-> %f:  %f\n",Q2,x,r,value,r*value/pow(1-R,2));
	value=value/pow(1-R,2);
	
        return(value);
//	return(r*value);
}

#if Z_INTEGRATE==1

double sigma_DIS(double x,double q2,double y, double * par) {

        double sigma_sum;

        double var[2];//={x,q2};
        *(var)=x;
        *(var+1)=q2;
        double *param[2];//={var,par};
        *(param)=var;
        *(param+1)=par;
        double res=0;
        //printf("x:%f , Q2: %f \t",x,q2);
        double ep=1.0e-6;
        simpson1dA(& sigma_r_integrand, param,0.0+ep, 0.9 ,100,&res);

        return res;
}
#elif ZINTEGRATE==0 
double sigma_DIS(double x,double q2,double y, double * par) {

	double sigma_sum;

	double var[2]={x,q2};
	double *param[2]={var,par};
	double res=0;

simpson2d(& sigma_integrand, param, 1.0e-10,1.0-1.0e-10,0.0,1.0,&res);
	
	return res;
}	
#endif
