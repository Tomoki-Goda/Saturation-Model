/////////////////////////////////////////////////////////////////////
// Gluon_GBW af
// af.init(par) //double* par sigma parameters 
// af(x,kt2,mu2) // 
//
// Dipole gluon
// Dipole_Gluon af
// af.init(par,integ) //Integrand integ; see below
// af(x,kt2,mu2)
//
//WW_Gluon
////////////////////////////////////////////////////////////////////
//
//GBW gluon
//
////////////////////////////////////////////////////////////////////
#ifndef DIPOLE_GLUON_HH
#define DIPOLE_GLUON_HH

#include<cmath>
#include<iostream>
#include<vector>
#include<string>
#include"control.h"
#include"control-default.h"
#include"constants.h"
#include"clenshaw.hh"
#include <gsl/gsl_sf.h>
#include"Levin.hh"
#include"gluon-integrand.hh"

typedef Gluon_Integrand INTEG ;

class Gluon_GBW{
	const double *sigpar;
	double sigma_0=0,lambda=0,x_0=0,mu02=0,thresh_power=0;
	double x=0;
	public:
		explicit Gluon_GBW(){
		}
		~Gluon_GBW(){}
		void init(const double *par);
		
	public:
		double operator()(const double  x,const double k2,double mu2);
};

/////////////////////////////////////////////////////////////////////
//
//Dipole gluon
//
/////////////////////////////////////////////////////////////////////
class Dipole_Gluon{
//class Dipole_Gluon{
		const double *par;
		INTEG *integrand;
		CCIntegral cc=CCprepare(64,"dipole",1,5);
	public: 
		Dipole_Gluon(INTEG & integ){
			integrand =& integ;
		}
		~Dipole_Gluon(){
			
		}
		void init(const double * const &par );
		void set_x(const double &x);
		double operator()(const double x,const double kt2,const double mu2);
		
		
};
#endif
