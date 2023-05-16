#include<wstp.h>
#include"gluons.hh"
Collinear_Gluon xgpdf;
//double Ag,lg;
double xg(const double x, const double QQ,const double A_g,const double l_g){
	//Collinear_Gluon xgpdf;
	return(xgpdf( x,  QQ, A_g,l_g));
}

Chebyshev1D_Collinear_Gluon xgpdf_cheb;
int xg_cheb_init(const double min, const double max,const double Ag, const double lg,double xx){
	static double x;
	x=xx;
	//Collinear_Gluon xgpdf;
	xgpdf_cheb.allocate(25);
	xgpdf_cheb.init(min,max, Ag,lg);
	xgpdf_cheb.set_x(x);
	return(0);
}

double xg_cheb(const double x, const double QQ){
	//Collinear_Gluon xgpdf;
	return(xgpdf_cheb( x,  QQ));
}
double alpha(double mu2 ){
	static double b0= ((double)(33 -2*NF))/(12*PI);
	return( 1/(b0* log(mu2/LQCD2)));//LQCD2 lambda_QCD ^2
}

int main(int argc,char** argv){
	return WSMain(argc, argv);
}


:Begin:
:Function: xg
:Pattern: XGluon[x_?NumberQ, q2_?NumberQ,ag_?NumberQ, lg_?NumberQ]
:Arguments: {N[x], N[q2],N[ag], N[lg]}
:ArgumentTypes: {Real, Real,Real, Real}
:ReturnType: Real
:End:

:Begin:
:Function: xg_cheb_init
:Pattern: XGluonChebInit[min_?NumberQ, max_?NumberQ,ag_?NumberQ, lg_?NumberQ,x_?NumberQ]
:Arguments: {N[min], N[max],N[ag], N[lg], N[x]}
:ArgumentTypes: {Real, Real,Real, Real,Real}
:ReturnType: Integer
:End:

:Begin:
:Function: xg_cheb
:Pattern: XGluonCheb[x_?NumberQ, q2_?NumberQ]
:Arguments: {N[x], N[q2]}
:ArgumentTypes: {Real, Real}
:ReturnType: Real
:End:

:Begin:
:Function: alpha
:Pattern: alpha[q2_?NumberQ]
:Arguments: {N[q2]}
:ArgumentTypes: { Real}
:ReturnType: Real
:End:
