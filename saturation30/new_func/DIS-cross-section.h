#include<iostream>

#include<cmath>
#include<string>
#include<vector>

#include"../main.h"

class cross_section{
	private:
		std::string type ="none";
		std::string model="GBW";
		double(*s_ptr)( double);
		double(*psi_ptr)(double,double);
		photon_psi psi_l('l');
		photon_psi psi_s('s');
		photon_psi psi_c('c');
		photon_psi psi_b('b');
		psi_L(double r,double x){
			return(  psi_l.psisq_f(r,z)+ psi_s.psisq_f(r,z));
		}	
		psi_LC(double r,double x){
                        return(  psi_l.psisq_f(r,z)+psi_s.psisq_f(r,z)+psi_c.psisq_f(r,z) );
                }
		psi_LCB(double r,double x){
                        return( psi_LC(r,z)+psi_b.psisq_f(r,z) );
                }

		psi_C(double r,double x){
                        return( psi_c.psisq_f(r,z) );
                }
		psi_B(double r,double x){
                        return( psi_b.psisq_f(r,z) );
                }
		
		
		double integrand(int *dim, double* x){
			return( x[0] *(*psi_prt)(x[0],x[1])*(*s_ptr)(x[0],x[1],Q2,par) );
		}
		double Q2,xmod,y;
		unsigned xbj_mod;
	public:
		cross_section(unsigned f, unsigned f_b, unsigned model_id, unsigned x_id ){
			set_type(f,f_b,model_id);
			xbj_mod=x_id;
		}
		void set_variables(double X,double Y,double q2,double* param){
                        /* Set the global variables */
                        //xmod = X;
			Q2=q2;
			par= param;
                        switch (xbj_mod){
                                case 0:
                                        xmod = X;
                                        break;
                                case 1:
                                        xmod = X*(1+4*m_fsq/Q2);
                                        break;
                        }
		}

		void set_type(unsigned f,unsigned f_b, int model_id){
			std::string* type_p=&type;
			std::string * model_p=&model;
			


			switch(model_id){
				case 0:
					*model_p="GBW";
					s_ptr=&sigma_gbw;
					break;
				case 1:
					*model_p="BGK";
					s_ptr=&sigma_bgk;
					break;
				case 2:
					*model_p="GBS";
					s_ptr=&sigma_gbs;
					break;
				}
			switch(f){
                                case 0:
                                        *type_p="Light";
					(psi_ptr)=&psi_L;
					//photon_psi psi_l('l');
			                //photon_psi psi_s('s');
                                        break;
                                case 1:
					if(f_b==0){
						*type_p="Light + Heavy";
						psi_ptr=&psi_LC;
						//photon_psi psi_l('l');
				                //photon_psi psi_s('s');
				                //photon_psi psi_c('c');
					}
					if(f_b==1){ 
						*type_p="Light + Heavy + B";
						psi_ptr=&psi_LCB;
						//photon_psi psi_l('l');
        			        	//photon_psi psi_s('s');
				                //photon_psi psi_c('c');
                				//photon_psi psi_b('b');
					}
					break;
                                case 2:
                                        *type_p="Charm";
					psi_ptr=&psi_C;
             				//photon_psi psi_c('c');
					break;
                                case 3:
                                        *type_p="Bottom";
					psi_ptr=&psi_B;
                			//photon_psi psi_b('b');
					break;
                                }
			
			}

		void show_parameters(std::string word){
			std::cout<< "Model: "<< *model_p<<std::endl;
			std::cout<< "Flavour: "<< *type_p<<std::endl;
			//std::cout<< "Something else"<<word<<std::endl;	
		}

		


		double evaluate() {
			//Options are, 
			//* Flavour
			//* Model

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

			B[0]=Rmax_simps2d;

			result=simps2d(A,B,eps_simps2d,&integrand);
		}
};

//////////////////////////////////////////////////////////////////////////////////////////////
//function
//////////////////////////////////////////////////////////////////////////////////////////////
/*cross_secton sigma_f(flavour, fl_beauty, model);

double sigma(double X, double Q, double Y, double *par){
	return( sigma_f.evaluate(X,Q,Y, par) );
}
*/


