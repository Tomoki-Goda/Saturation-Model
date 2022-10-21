//------------------------------------------------------------------------------
//
// Copyright (c) 2012, Sebastian Sapeta
//
//  Modified for GKS gluon by Tomoki Goda, 21/10/2022.
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
#ifndef __INTERPOLATION_HH__
#define __INTERPOLATION_HH__
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

using namespace std;

//------------------------------------------------------------------------------
/// \class BaseInterpolation
/// Base class for interpolation
//------------------------------------------------------------------------------
class BaseInterpolation {
  public:
  /// default constructor
  BaseInterpolation() {}
  /// default destructor
  virtual ~BaseInterpolation() {}

  /// grid limits in form of vector of pairs: [[x1min,x1max], [x2min,x2max],..]
  virtual vector<pair <double, double> > grid_limits() const { 
    return vector<pair <double, double> >(); }

  /// return the interpolated value corresponding to the vector (x1,x2,x3,..)
  virtual double interp (const vector<double> &x) { return -1.0; }
  
  /// return status of the last evaluation
  ///
  /// It needs to be implemented in each derived class, the idea is that the
  /// variable is recomputed after each evaluation. 
  /// Its binary format is such that the
  /// rightmost bit corresponds to the first variable from the vector x of the
  /// interp(const vector<double> &x) so for example: 11 (decimal 3) the first
  /// two elements were out of range, 10 (decimal 2) the second element was out
  /// of range, 0 (decimal 0) none element out of range (in other words
  /// success).
  unsigned int eval_status() {return _eval_status;}

  protected:
  unsigned int _eval_status; ///< store info about status of the last evaluation

};


//------------------------------------------------------------------------------
/**
 \class Interpolation2D
 2-dim interpolation - linear or cubic spline 


 The class uses 1d gsl interpolation function.

 The evalaution is done with the function interp(const vector<double> &x) and
 the status of the last evaluation can be checked with eval_status().

 Example of usage:

 for (unsigned i=0;i<n;i++) {         
   for (unsigned j=0;j<m;j++) {      
       ...
       vec_x1[i] = arg1;
       vec_x2[j] = arg2;
       matrix_A[i][j] = val;
   }
 }

 Interpolation2D int2d(vec_x1, vec_x2,  matrix_A, *gsl_interp_cspline);
 ...
 int2d.interp(vec_x);
 \endverbatim 
**/ 
class Interpolation2D : public BaseInterpolation {

  public:
  /// constructor -
  /// read the input and interpolate in the innermost dimension x2
  ///
  /// \param x1vv, x2vv   input vectors with argument data; entries in both
  ///                     vectors are required to grow monotonically
  /// \param ym           input array with values data
  /// \param int_type     type of interpolation, it takes the values
  ///                     values *gsl_interp_linear or *gsl_interp_cspline
  Interpolation2D(const vector<double> &x1vv, const vector<double> &x2vv,
                 double *ym, 
		 gsl_interp_type int_type = *gsl_interp_linear);

  /// destructor		
  ~Interpolation2D();
  
  /// return grid limits in form of 2d vector of pairs: 
  /// [[x1min,x1max], [x2min, x2max]]
  /// 
  /// Min and max values in each dimension can be obtained by using
  /// first  or second pair data members.
  vector<pair <double, double> > grid_limits() const;

  /// return value from 2d spline interpolation at (x1,x2)
  double interp (const vector<double> &x);

  private:
  // type of interpolation
  gsl_interp_type _int_type; 
  // sizes of the multi_array in each of two dimensions (x1, x2)
  unsigned int n1, n2;
  double* y;
  vector<gsl_interp_accel*> acc;
  vector<gsl_spline*>  spline;
  // vectors storing x1 and x2 grids
  vector<double> x1v,x2v;
};

//------------------------------------------------------------------------------
//
/// \class CubicSpline3D
///  3-dim interpolation with cubic spline 
///
/// Interpolation of 3D function with the cubic spline method.
///
/// @todo
///   - either document this class properly or think about excluding it from the
///   distribution
//------------------------------------------------------------------------------
class CubicSpline3D : public BaseInterpolation {

  public:
  /// constructor
  /// get all input and interpolate in the innermost dimension x3
  CubicSpline3D(const vector<double> &x1vv, const vector<double> &x2vv, 
                const vector<double> &x3vv, double *ym);
		
  ~CubicSpline3D();

  //double interp (double x3) {
  //      int i = 0, j = 0;
  //      return gsl_spline_eval (spline[i][j], x3, acc[i][j]);
  //}
  
  /// return value from 3d spline interpolation
  /// get result for interpolation in x3 and use it to do two consecutive
  /// interpolations in x2 and x1

/// Min and max values in each dimension can be obtained by using
  /// first  or second pair data members.
  vector<pair <double, double> > grid_limits() const;
  double interp (const vector<double> &x);
  

  private:
  // sizes of the multi_array in each of three dimensions
  unsigned int n1, n2, n3;
  double* y;
  vector<vector<gsl_interp_accel*> >acc;
  vector<vector<gsl_spline*> > spline; 
  // vectors storing x1, x2 and x3 grids
  vector<double> x1v,x2v,x3v;
};


#endif // __INTERPOLATION_HH__
