#ifndef DCLENSHAW2_HH
#define DCLENSHAW2_HH 1
#include<math.h>
#include<stdio.h>
#include<string>
#include"./Kahn.hh"
//#include<stdlib.h>

#ifndef PI
#define PI 3.141592653589793238462643383279502884197
#endif
//
//To use this, one needs to define cc_integrand in a drived class. 
// then cc_integrate integrates cc_integrate between a and b with extra parameter par 
// to the accuracy eps and absolute accuracy Aeps
//

class Clenshaw_Curtis{
	private:
		std::string tag;
		int N=128;
		double wfull[129], whalf[65], x[129];
		int max_rec=5;
		int InitDiv=1;
		inline int t_pos(int i,int j, int n){
			return(n/4-abs(n/2-i*j+((i*j)/n)*n));
		}
		inline int sign(int i){
			return( (i==0)?(1): (i/abs(i)));
		}

		double fixed_cc(const void* par,const double smin,const double smax, Kahn &full, Kahn &half)const;
	protected:		
		int cc_init(const int N);
		int cc_init(const int N, const std::string& name,int max_rec);
		int cc_init(const int N, const std::string& name,int InitDiv,int max_rec);
	public:
		explicit Clenshaw_Curtis(){
		}
		~Clenshaw_Curtis(){
		}
		virtual double cc_integrand(const double,const void*)const=0;
		double cc_integrate(const void* par,const double a, const double b, const double eps, const double Aeps)const;
};
#endif // DCLENSHAW2_HH
