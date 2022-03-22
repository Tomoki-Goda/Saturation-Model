#include<iostream>
#include<fstream>
#include<cmath>

#include <string>
#include "gluons.h"
#include "float.h"
#include "chebyshev.h"
#include "chebyshev3.h"
#include "/home/tomoki/Tools/Cuba-4.2/cuba.h"

/********************************************************************************/
/*                               Photon Wave Function                           */
/********************************************************************************/

class PhotnWF {
/*photon wave function */
/*See for example J. R. Forshaw, et al 1999*/
	private:
		
		
	public:
		/* light, strange ,etc  number coded*/
		int particle_type;
		
		/* parameters & control */
		double	m_fsq=0;
		double	Qsq=1; 
		double	y ;
		int	dataform = 0;
		float	nf = 5.0/6.0; //seemingly sum of fraction charge squared times 3/2?
		float	norm=1;   //in main.h  ~0.8 
		
		
		
		
		double psi2(double r ,double z){
		/*Photon wave function |psi|^2 */
		
		  	double     z_bar =  pow(z,2) + pow(1-z,2) ;
		  	double     y_bar =  pow( y ,2) / (1 + pow(1-y,2));
		  	double   Qsq_bar =  z * (1-z) * Qsq + m_fsq;
		  	double      Qsq2 =  sqrt(Qsq_bar) * r;
		  	double bessel_k0 = dbesk0_(&Qsq2);
		  	double bessel_k1 = dbesk1_(&Qsq2);
			double value;
   	
		
			switch(dataform) {
				case 0: /* Reduced cross section form */  
					value =  norm * nf * Q2 * (z_bar * Qsq_bar* pow( bessel_k1 , 2 ) + ( m_fsq + (1-y_bar) * 4 * Q2 * pow( z , 2 ) * pow( 1-z , 2) ) * pow(bessel_k0,2) ) ;
				break;
				case 1: /* F_2 form */  
					value =  norm * nf * Q2 * (z_bar * Qsq_bar* pow( bessel_k1 , 2 ) + ( m_fsq + 4 * Q2 * pow( z , 2 ) * pow( 1-z , 2 ) ) *  pow( bessel_k0 ,2) );
					//value =  norm * nl * Q2;	
				break;
				case 2: /* F_L form */
					value =  norm * nf * Q2 * 4 * Q2 * pow( z , 2 )* pow( 1-z , 2 ) * pow( bessel_k0 , 2 );
				break;
			}
			return(value);
		};


};

/*****************************************************************************************************************/

/*****************************************************************************************************************/

class CrossSection {
	private:
	
	public:
		
};


/****************************************************************************************************************/
/*                                      Dipole crosssection                                                     */
/****************************************************************************************************************/



/*
class DipoleCorossSection{
	private:
		
		double exparg( )
		
	public:
		
		double	sigma0	= 1;
		double	exparg	= 1;
		int	model	= 0;
		
		switch(model){
			case 0:
				
		} 
				
		
		double SigmaDP(double r, double x){
		double out = sigma0 * (1- exp( exparg ));
		
		return(out);
		
		};
		
};
*/

/*******************************************************************************************************/
/*                                    miscelaneous functions etc                                       */
/*******************************************************************************************************/


/**************************************************************/
/*        The strong running coupling alpha_s                 */
/**************************************************************/
double alpha_s (double Q) {

    double b0;

    b0 = (33-2*n_f)/(12*pi);
    
    return 1/(b0*log(Q/Lambda2));  /* Q means Q^2 [GeV^2]*/

}


