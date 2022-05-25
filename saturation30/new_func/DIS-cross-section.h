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
	}
	return (x * (1 + 4 * m_fsq / Q2));
}

double sigma_integrand(double R,double z,double * par[2]){

	//R is r/(1+r) such that integration between (0, inf) -> (0,1).
	double r=R/(1-R);
	double jacob=pow(1-R,-2);	
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

	return(jacob*r*value);
}


double sigma_r_integrand(double R,double * par[2]){

        //R is r/(1+r) such that integration between (0, inf) -> (0,1).
        double r=R/(1-R);
        double jacob=pow(1-R,-2);
        double x=(*(*par));
        double Q2=(*( (*par)+1 ));

        double* param=( *(par+1) );
        double value=0;
        //double x_mod[4];

        #if FLAVOUR==0  
                value=SIGMA(r,mod_x(x,Q2,'l'), Q2, param) * psisq_z_int(r, Q2, 'l');
                value+=SIGMA(r,mod_x(x,Q2,'s'), Q2, param) * psisq_z_int(r, Q2, 's');
        #elif FLAVOUR==1
                        value=SIGMA(r,mod_x(x,Q2,'l'), Q2, param) * psisq_z_int(r, Q2, 'l');
                        value+=SIGMA(r,mod_x(x,Q2,'s'), Q2, param) * psisq_z_int(r, Q2, 's');
                        value+=SIGMA(r,mod_x(x,Q2,'c'), Q2, param) * psisq_z_int(r, Q2, 'c');
        #elif FLAVOUR==2
                        value=SIGMA(r,mod_x(x,Q2,'l'),Q2, param) * psisq_z_int(r, Q2, 'l');
                        value+=SIGMA(r,mod_x(x,Q2,'s'), Q2, param) * psisq_z_int(r, Q2, 's');
                        value+=SIGMA(r,mod_x(x,Q2,'c'), Q2, param) * psisq_z_int(r, Q2, 'c');
                        value+=SIGMA(r,mod_x(x,Q2,'b'), Q2, param) * psisq_z_int(r, Q2, 'b');
        #endif

        return(jacob*r*value);
}

#if Z_INTEGRATE==1

double sigma_DIS(double x,double q2,double y, double * par) {

        double sigma_sum;

        double var[2]={x,q2};
        double *param[2]={var,par};
        double res=0;

simpson1d(& sigma_r_integrand, param,1.0e-10,1.0-1.0e-10,&res);

        return res/PI;
}
#elif ZINTEGRATE==0 
double sigma_DIS(double x,double q2,double y, double * par) {

	double sigma_sum;

	double var[2]={x,q2};
	double *param[2]={var,par};
	double res=0;

simpson2d(& sigma_integrand, param, 1.0e-10,1.0-1.0e-10,0.0,1.0,&res);
	
	return res/PI;
}	
#endif
