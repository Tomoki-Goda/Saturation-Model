#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>

#include"../main.h"

class photon_psi{
	//photon wave function.
	// use set_type(char 'f') to set which type of 'l', 's', 'c', 'b'.
	// 
	// main function is psisq(double r, double z),
	// with Global variable Q2 used.
	//
	//show_parameters() to show mass and charge_sum used.
	//
	private:
		std::string type= "light";
		double mass=m_fsq;
                double charge_sum=5.0/6.0;
	public:
		//constructor, default is light flavours.
		photon_psi(){
		set_type('l');
		}
		
		void set_type(char flavourtype){
			//pointers to mass and charge_sum, global in the class.
			double* mass_p=&mass;
			double* charge_p=&charge_sum;
			std::string* type_p=&type;
			//mass are defined in main.h
			switch(flavourtype){
				case 'l':
					*charge_p=5.0/6.0;
					*mass_p=m_fsq;
					*type_p="Light";
					break;
				case 's':
					*charge_p=1.0/6.0;
					*mass_p=m_s;
					*type_p="Strange";
					break;
				case 'c':
					*charge_p=2.0/3.0;
                                        *mass_p=m_ch;
					*type_p="Charm";
                                        break;
				case 'b':
					*charge_p=1.0/6.0;
                                        *mass_p=m_b;
					*type_p="Bottom";
                                        break;
			}
		}
		
		 void show_parameters(){
			std::cout<<"Type :" <<type<<std::endl;
                        std::cout<<"mass :"<<mass<<std::endl;
                        std::cout<<"charge sum :"<<charge_sum<<std::endl<<std::endl;
                }
		
		//MAIN FUNCTION.
		//
		double psisq (double r, double z)  {
			//!! Q2 is somehow global variable...
			double	z_bar =  z*z+(1-z)*(1-z);
			//double     y_bar =  (y*y)/(1+(1-y)*(1-y));
			double	Qsq_bar =  z*(1-z)*Q2+mass;
			double	Qsq2 =  sqrt(Qsq_bar)*r;
			double	norm_f;
			//bessels functions, declared in main.h
			double	bessel_k0 = dbesk0_(&Qsq2);
			double	bessel_k1 = dbesk1_(&Qsq2);
			
			double	value;
	
			
			switch (dataform) {
				case 1: // F_2 form
					value = (charge_sum)*norm*Q2*(z_bar*Qsq_bar*bessel_k1*bessel_k1
						    + ( mass +4*Q2*z*z*(1-z)*(1-z))*bessel_k0*bessel_k0);
					  
				break;
				case 2: // F_L form 
					value = (charge_sum)*norm * Q2*4*Q2*z*z*(1-z)*(1-z)*bessel_k0*bessel_k0;
		        		break;
		    	}
		
			if (photo ==1) {
				norm_f = (charge_sum)*(charge_sum);
				value = norm_f*z_bar*Qsq_bar*bessel_k1*bessel_k1;
			}
			
			return (value);
		}

	
};







