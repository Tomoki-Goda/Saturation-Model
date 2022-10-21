//------------------------------------------------------------------------------
//
// Copyright (c) 2012, Sebastian Sapeta
//  Modified for GKS gluon by Tomoki Goda, 21/10/2022.
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
#include "UnintegratedGluonGKS.hh"
#include "InterpolationGKS.hh"

using namespace std;

int main(int argc, char ** argv) {


  /// Set interpolation type: 
  ///   *gsl_interp_linear 
  ///   *gsl_interp_cspline
  gsl_interp_type int_type =  *gsl_interp_cspline;

  /// Choose the grid:
  ///    grids/KSlinear.dat
  ///    grids/KSnonlinear.dat
  ///    grids/KShardscalelinear.dat
  ///    grids/KShardscalenonlinear.dat
  ///    grids/KSsudFqgnonlinear.dat
  ///    grids/KSsudFggnonlinear.dat
  
  //string gridfile =  "grids/KSlinear.dat";
  //string gridfile =  "grids/KSsudFggnonlinear.dat";
  //string gridfile =  "grids/KShardscalelinear.dat";
  string gridfile =  "grids/GBWS-gluon.dat";
  //printf("Current selection= %s\n",gridfile);
  string response;
  bool res=false;
  while(!res){
        std::cout<<"Current selection="<<gridfile<<std::endl;
	std::cout<<"This File ? type file name or type Y to use the currently chosen file:\n";
  	//printf("This File ? type file name or type Y to use the currently chosen file:\n" );
	//scanf("%s",respose);
	std::cin>>response;
	if((response.compare("y")* response.compare("Y"))==0){
		res=true;
	}else{
		gridfile=response;
	}
  }
  std::cout<<"Start computing"<<std::endl;	  

  /// Created updf object
  UnintegratedGluon updf = UnintegratedGluon(gridfile,int_type);


  /// Use interpolation
  double kt = 1.0;
  double mu =  1.0;
  cout << "# --------------------------------------------------------" << endl;
  cout << "# "<< gridfile           << endl;
  cout << "# mu = " << mu << " GeV" << endl;
  cout << "# kt = " << kt << " GeV" << endl;
  cout << "# --------------------------------------------------------" << endl;
  cout << "#   x       xg(x,kt)"     << endl; 
  cout << "#" << endl;
  //cout << "# --------------------------------------------------------" << endl;
  unsigned int npoints = 20;
  for (unsigned i = 0; i<npoints; i++) {
      double logx = -15.0 + 13.0*i/npoints;
      // 3rd argument is relevant only for the KShardscale grids
      cout << exp(logx) << "  " << updf.xg(logx,log(kt*kt), log(mu*mu)) << endl;
  }

}
