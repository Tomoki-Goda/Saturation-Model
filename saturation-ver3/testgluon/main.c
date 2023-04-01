#include<stdio.h>
#include<math.h>
#include"./gluons.h"

//extern int simps_(double*,double*,double*,double*,double*,
  //             double(*)(double*),double*,double*,double*,double*);
int main(){
	set_xg_parameter(1,0.5995);
	double val=xgpdf(1.0e-8,6.0e+3);
	printf("%.3e\n",val);
	return 0;
}
