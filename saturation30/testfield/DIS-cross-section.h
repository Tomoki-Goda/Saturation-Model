#include<iostream>
#include<cmath>
#include<string>
#include<vector>

#include"../main.h"

class cross_section{
	private:
		std::string type ="none";
		std::string model="GBW";	
	public:
		cross_section(){
			set_type('l',0);
		}
		void set_type(char f, int model_id){
			std::string* type_p=&type;
			std::string * model_p=&model;
			double * integrand_p;
					
			switch(model_id){
				case 0:
					*model_p="GBW";
					integrand_p=&()

					break;
				case 1:
					*model_p="BGK";
					break;
				case 2:
					*model_p="GBS";
					break;
				}
			switch(f){
                                case 'l':
                                        *type_p="Light";
                                        break;
                                case 's':
                                        *type_p="Strange";
                                        break;
                                case 'c':
                                        *type_p="Charm";
                                        break;
                                case 'b':
                                        *type_p="Bottom";
                                        break;
                                }

			}
			


		
		void show_parameters(std::string word){
			std::cout<< "Model: "<<model<<std::endl;
			std::cout<< "Flavour: "<<type<<std::endl;
			std::cout<< "Something else"<<word<<std::endl;	
		}




		double sigma_f(double X, double Q, double Y, double *par) {
			
			/* Integration limits */
			double A[2] = {amin,0.0};/* Lower limits of integration  {r,z} */
			double B[2] = {Rmax,0.5};       /* Upper limits of integration  {r,z} */
			
			
    			/* dadmul integration variables */
			int    dim = 2;            /* Dimension of integral */
			int    minpts = 1.7e+01;
			int    maxpts = 1.0e+05;
			//int    maxpts = 1.0e+04;
			double eps = sigmaEPS;     /* Relative accuracy */
			//int    iwk = 2100;
			int    iwk = 21000;
			double wk[iwk];
			double result;             /* The result of integration */
			double relerr;             /* The relative accuracy of result */
			int    nfnevl;             /* The number of function eval. performed */
			int    ifail;              /* If '0' normal exit relerr < eps */
			//double temp;

			/* Set the global variabels according to sigma() arguments */
			//xmod = X;
			switch (xbj_mod){
				case 0:
					xmod = X;
					break;
				case 1:
					xmod = X*(1+4*m_fsq/Q);
					break;
			}

			Q2  = Q;
			y   = Y;

		
			B[0]=Rmax_simps2d;
			result=simps2d(A,B,eps_simps2d,*integrand_function);
		}
};
