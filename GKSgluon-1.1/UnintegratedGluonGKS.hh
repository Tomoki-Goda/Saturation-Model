//------------------------------------------------------------------------------
//
// Copyright (c) 2012, Sebastian Sapeta
//
//  This is part of KSgluon package.
//  Modified for GKS gluon by Tomoki Goda, 21/10/2022.
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
#ifndef __UNINTEGRATEDGLUON_HH__
#define __UNINTEGRATEDGLUON_HH__

#include <string>
#include <vector>
#include "InterpolationGKS.hh"
using namespace std;

//------------------------------------------------------------------------------
/// enums for gluon names
enum UnintegratedGluonName {
  grid2d,  // gluon interpolated from 2d grid
  grid3d   // gluon interpolated from 3d grid
};


//------------------------------------------------------------------------------
/**
\class UnintegratedGluon

Computing unintegrated gluon by interpolating an input grid

The form of the 2d grid is assumed to be:
\verbatim
 # logx   logkt2    xg
    i1     j1      val1
    ...    ...     ...
    i1     jn      valn
    i2     j1      valn+1
    ...    ...     ...
    i2     jn      valn*m     
\endverbatim

and both variables have to grow monotonically.

IMPORTANT: normalization of gluon

In case of problems in gluon evaluation a corresponding error is written to the
\p unintegratedgluon.err file which is always produced by this class. The
possible error codes are

 - Err001 : Evaluation outside the domain of the grid.
 - Err002 : Negative value of the gluon.

In each of the above cases the return value of the xg(double logx, double logkt2) function is set to zero.

@todo 
 - add description of normalization, perhaps remove functions _norm_f/F

**/ 

class UnintegratedGluon {
  public:
  /// base constructor taking the grid at input
  /// \param filename name of the input grid file
  /// \param int_type type of interpolation; it takes gsl
  /// defined values *gsl_interp_linear or *gsl_interp_cspline
  UnintegratedGluon (string filename, gsl_interp_type int_type);
  /// return size of the grid in each dimension
  //pair<unsigned int,unsigned int> grid_size() const;
  vector<int> grid_size() const;
  /// return grid limits in form of 2d vector of pairs: 
  /// [[x1min,x1max], [x2min, x2max]]
  /// 
  /// Min and max values in each dimension can be obtained by using
  /// first  or second pair data members.
  vector<pair <double, double> > grid_limits() const {
    return _updf->grid_limits();}
  /// value of the unintegrated gluon
  /// \param logx    ln(x) where x is the longitudinal momentum fraction of the
  ///                gluon
  /// \param logkt2  ln(kt^2) where kt is the transverse momentum of the gluon
  /// \param logmu2  ln(\mu^2) where \mu is the factorization scale; 
  ///                this parameter is relevant only for some gluons like 
  ///                e.g. the KShardscale gluon; it also implies that 3d grid 
  ///                is used as input
  double xg(double logx, double logkt2, double logmu2 = 0.0);
  virtual ~UnintegratedGluon ();
  private:
  // function determining the size of the 2d grid
  void _get_2dgrid_dimensions(string filename);
  // function determining the size of the 3d grid
  void _get_3dgrid_dimensions(string filename);
  /// set updf as a function interpolated from 2d grid
  ///
  /// the form of the grid is assumed to be:
  ///     logx   logkt2   xg
  ///     i1     j1
  ///     ...    ...
  ///     i1     jn
  ///     i2     j1
  ///     ...    ... 
  ///     i2     jn              
  void _set_from_2dgrid(string filename);
  void _set_from_3dgrid(string filename);
  // functions setting gbw pdf
  void _set_gbw();
  UnintegratedGluonName _determine_grid_type(string filename);
  unsigned int _n1, _n2, _n3;
  ofstream *_errorfile;
  static double _norm_F(const vector<double> &x);
  BaseInterpolation *_updf;
  UnintegratedGluonName _name;
  // pointer to the normalization function for pdf
  double (*_norm)(const vector<double> &x);
  gsl_interp_type _int_type;
  double _some_eps;
};


#endif // __UNINTEGRATEDGLUON_HH__
