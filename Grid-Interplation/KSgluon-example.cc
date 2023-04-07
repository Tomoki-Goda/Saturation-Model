//------------------------------------------------------------------------------
//
// Copyright (c) 2012, Sebastian Sapeta
//
//  This is part of KSgluon package.
//
//  KSgluon is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  KSgluon is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet; if not, write to the Free Software
//  Foundation, Inc.:
//      59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//------------------------------------------------------------------------------

#include <iostream>
#include <string>
#include "UnintegratedGluonKS.hh"
#include "InterpolationKS.hh"

using namespace std;



//UnintegratedGluon prepare_gluon(std::string gridfile,int dimension) {
int main(int argc ,char** argv){ 
	gsl_interp_type int_type =  *gsl_interp_cspline;
	std::string gridfile=argv[1];
	UnintegratedGluon updf = UnintegratedGluon(gridfile,atoi(argv[2]),int_type);
	FILE* file=fopen(argv[3],"w");
	
	double x=1.0e-3;
	double mu2=17*17;
	double k2=0;
	for(int i=0;i<=100;++i){
		k2=exp(8*((double)i)/99);
		fprintf(file, "%.5e\t%.5e\n",k2,updf.xg(log(x),log(k2),log(mu2)));
	}
	fclose(file);

	return 0;
}
