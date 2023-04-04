#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<cmath>
#include"./Functions/clenshaw.hh"
//#define PI 3.14159265358979323846264338327950288419716939937510

int k=100;

class Integrand{
	public:
		/*double operator()(double r,void*){
			return(r*exp(-r));
		}*/
		double operator()(double r,void*p){
			double *i=(double*)p;
			return(i[1]*cos( r * i[0])+pow(i[1],2)*sin(r*i[0]*2));
			//return( std::cyl_bessel_j(1,k*r) );
		}
			
};

int main (){
	Integrand integrand;
	CCIntegral cc=CCprepare(64,"name",1,5);
	double par[]={50,1000};
	double res;

	res =dclenshaw<Integrand,void*>(cc,integrand,(void*)par,0,2*PI,1.0e-14,1.0e-15);
	printf("1:  %.3e\n",res);


	double min=0 ,max=2*PI;
	double sum=0;
	double scale=5*(max-min)/par[0];
	for (int i=0;i<par[0];i++){	
		max+=scale;

		res =dclenshaw<Integrand,void*>(cc,integrand,(void*)par,min,max,1.0e-14,1.0e-15);
		printf("2.%d:%.3e \n",i,res);
		min=max;
		sum+=res;
	}

	printf("2: %.3e\n",res);
	
	return 0;
}
