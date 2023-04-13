#include<iostream>
#include"gluons.hh"

extern "C" double xgpdf_(const double* x, const double* QQ,const double* A_g, const double* lambda_g );


int main(){
	Collinear_Gluon xg;
	Chebyshev1D_Collinear_Gluon cxg(30);
	double x,mu2;
	double val1,val2,val3;
	const double Ag=1.1,lg=0.08;
	cxg.init(1,2.0e+8,Ag,lg);
	for(int j=0;j<50; ++j){
		x=1.0e-6*pow(1.0e+6,((double)j)/49);
		cxg.set_x(x);
	for(int i=0;i<50;++i){	
		mu2=2*pow(1.0e+8/2,((double)i)/49);
		val1=xgpdf_(&x,&mu2,&Ag,&lg);
		val2=xg(x,mu2,Ag,lg);
		val3=cxg(x,mu2,Ag,lg);
		if(fabs((val1-val2)/(val1+val2))>1.0e-6||fabs((val1-val3)/(val1+val3))>1.0e-6){
			printf("FORTRAN=%.2e\tdiff=%.2e,\t%.2e\n",val1,val1-val2,val1-val3);
		}
	}
	}
		
		
	return 0;
}
