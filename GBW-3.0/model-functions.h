#include<iostream>
#include<fstream>
#include<cmath>

#include <string>
#include "gluons.h"
#include "float.h"
#include "chebyshev.h"
#include "chebyshev3.h"
#include "/home/tomoki/Tools/Cuba-4.2/cuba.h"



class photnWF {
/*photon wave function */
	private:
		
		
	public:
		/* light, strange ,etc  number coded*/
		int particle_type;
		
		/* parameters & control */
		double	m_fsq=0;
		double	Q2=1; 
		double	y ;
		int	dataform = 0;
		float	nf = 5.0/6.0;
		float	norm=1;
		
		
		
		
		double psi2(double r ,double z){
		/*Photon wave function |psi|^2 */
		
		  	double     z_bar =  z * z + (1-z) * (1-z);
		  	double     y_bar =  (y * y) / (1 + (1-y) * (1-y));
		  	double   Qsq_bar =  z * (1-z) * Q2 + m_fsq;
		  	double      Qsq2 =  sqrt(Qsq_bar) * r;
		  	double bessel_k0 = dbesk0_(&Qsq2);
		  	double bessel_k1 = dbesk1_(&Qsq2);
			double value;
   	
		
			switch(dataform) {
				case 0: /* Reduced cross section form */  
					value =  norm * nf * Q2 * (z_bar*Qsq_bar*bessel_k1 * bessel_k1 + (m_fsq+(1-y_bar)*4*Q2*z*z*(1-z)*(1-z)) * bessel_k0 * bessel_k0);
				break;
				case 1: /* F_2 form */  
					value =  norm * nf * Q2 * (z_bar*Qsq_bar*bessel_k1 * bessel_k1 + (m_fsq+4*Q2*z*z*(1-z)*(1-z)) * bessel_k0 * bessel_k0);
					//value =  norm * nl * Q2;	
				break;
				case 2: /* F_L form */
					value =  norm * nf * Q2*4*Q2*z*z*(1-z)*(1-z)*bessel_k0*bessel_k0;
				break;
			}
			return(value);
		};


};

/*****************************************************************************************************************/

/*****************************************************************************************************************/

class crossSection {
	private:
	
	public:
		
};


/****************************************************************************************************************/
/*                                      Dipole crosssection                                                     */
/****************************************************************************************************************/


/*******************************************************************************
* The dipole-proton GBW cross section 
*******************************************************************************/
double sigma_gbw (double r) {

   return sigma_0*(1-exp(-0.25*r*r*pow(x_0/xmod,lambda)));
}


/***********************************************************************************************************/
/*                                   The strong running coupling alpha_s                                   */
/***********************************************************************************************************/
double alpha_s (double Q) {

    double b0;

    b0 = (33-2*n_f)/(12*pi);
    
    return 1/(b0*log(Q/Lambda2));  /* Q means Q^2 [GeV^2]*/

}


