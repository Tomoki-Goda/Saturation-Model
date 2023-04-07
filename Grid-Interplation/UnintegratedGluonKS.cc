#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include "UnintegratedGluonKS.hh"
using namespace std;


//-----------------------------------------------------------------------------
/*
UnintegratedGluon::UnintegratedGluon(string filename, gsl_interp_type int_type)
				     : _int_type(int_type), _some_eps(1e-10) {

  _name = _determine_grid_type(filename);

  _errorfile = new ofstream("unintegratedgluon.err");

  // determine grid dimensions (in logx and logkt2 directions) 
  // and set gluon to the interpolation of the 2d grid
  switch (_name) {
    case grid2d:
      _get_2dgrid_dimensions(filename);
      _set_from_2dgrid(filename);
      break;
    case grid3d:
      _get_3dgrid_dimensions(filename);
      _set_from_3dgrid(filename);
    default:
      break;
  }

}*/
UnintegratedGluon::UnintegratedGluon(string filename,int dimension, gsl_interp_type int_type)
				     : _int_type(int_type), _some_eps(1e-10) {

 // _name = _determine_grid_type(filename);

  _errorfile = new ofstream("unintegratedgluon.err");

  // determine grid dimensions (in logx and logkt2 directions) 
  // and set gluon to the interpolation of the 2d grid
  switch (dimension) {
    case 2:
      _get_2dgrid_dimensions(filename);
      _set_from_2dgrid(filename);
      break;
    case 3:
      _get_3dgrid_dimensions(filename);
      _set_from_3dgrid(filename);
      break;
    default:
    	printf(" Dimension has to be 2 or 3. dimension = %d\n",dimension );
      break;
  }

}


//-----------------------------------------------------------------------------
UnintegratedGluon::~UnintegratedGluon () {
  delete _updf; 
  _errorfile->close();
  delete _errorfile;
}

//-----------------------------------------------------------------------------
// function determining the size of the 2d grid
void UnintegratedGluon::_get_2dgrid_dimensions(string filename) {

  ifstream indata;
  indata.open(filename.c_str());
  if(!indata) { cerr << "Error: updf grid file could not be opened" << endl;
     exit(1); }

  // read all lines in the file; search for the turn of the first
  // variable, this gives n1, then continue reading to the the total
  // number of lines
  double x1, x2, f;
  int nlines = 0;
  bool count2 = true;
  vector<double> x1_vec;
  while(indata >> x1 >> x2 >> f) {
    ++nlines;
    if (count2) {
      x1_vec.push_back(x1);
      if (x1_vec.size()>1) {
        if(abs(x1_vec[x1_vec.size()-1]- x1_vec[x1_vec.size()-2])>_some_eps) {
	  count2 = false;
	  continue;
        }
      }
    }
  }

  // set the grid size in the argument variables, check if the grid is
  // a square in 2d  and print its size
  _n2 = x1_vec.size()-1;
  _n1 = nlines/ _n2;
  if (nlines % _n2) {
     cerr << "Error: updf grid is not a square 2d grid! (nlines/n1 = ";
     cerr << nlines << "/" << double(_n2) << " = ";
     cerr << double(nlines)/double(_n2) << ")"<< endl;
     exit(1);
  } 
  cout << "# grid2d size (" << _n1 << ", " << _n2 << ")" << endl;
  indata.close();
}

//-----------------------------------------------------------------------------
// function determining the size of the 3d grid
void UnintegratedGluon::_get_3dgrid_dimensions(string filename) {

  ifstream indata;
  indata.open(filename.c_str());
  if(!indata) { cerr << "Error: updf grid file could not be opened" << endl;
     exit(1); }

  // read all lines in the file; search for the turn of the first
  // variable, this gives n1, then continue reading to the the total
  // number of lines
  double x1, x2, x3, f;
  int nlines = 0;
  bool count2 = true;
  bool count3 = true;
  vector<double> x1_vec, x2_vec;
  while(indata >> x1 >> x2 >> x3 >> f) {
    ++nlines;
    if (count2) {
      x1_vec.push_back(x1);
      if (x1_vec.size()>1) {
        if(abs(x1_vec[x1_vec.size()-1]- x1_vec[x1_vec.size()-2])>_some_eps) {
          count2 = false;
          continue;
        }
      }
    }
    if (count3) {
      x2_vec.push_back(x2);
      if (x2_vec.size()>1) {
        if(abs(x2_vec[x2_vec.size()-1]- x2_vec[x2_vec.size()-2])>_some_eps) {
          count3 = false;
          continue;
        }
      }
    }
  }

  //cout << x1_vec.size()-1 << " " << x2_vec.size()-1 << endl;

  // set the grid size in the argument variables, check if the grid is
  // a square in 2d  and print its size
  _n3 = x2_vec.size()-1;
  _n2 = (x1_vec.size()-1)/_n3;
  _n1 = nlines/ (_n2*_n3);
  if (nlines % (_n2*_n3)) {
     cerr << "Error: updf grid is not a square 2d grid! (nlines/n1 = ";
     cerr << nlines << "/" << double(_n2) << " = ";
     cerr << double(nlines)/double(_n2) << ")"<< endl;
     exit(1);
  } 
  cout << "# grid3d size (" << _n1 << ", " << _n2 << ", " << _n3 << ")" << endl;
  indata.close();
}


//------------------------------------------------------------------------------
// return size of the grid in each dimension
vector<int> UnintegratedGluon::grid_size() const {
  vector<int> size_vec;
  size_vec.push_back(_n1);
  size_vec.push_back(_n2);
  size_vec.push_back(_n3);
  return size_vec;
}

//------------------------------------------------------------------------------
// set gluon from 2d grid as updf
void UnintegratedGluon::_set_from_2dgrid(string filename) {

  // read the grid input file
  ifstream indata;
  indata.open(filename.c_str());
  if(!indata) { cerr << "Error: updf grid file could not be opened" << endl;
     exit(1);}

  // vectors for passing argument values
  vector<double> vec_logx(_n1);
  vector<double> vec_logkt2(_n2);
  vector<double> vec_xg(_n2);

  // 2d array for passing function values at grid points
  double *matrix_xg = new double[_n1*_n2];

  // Read the grid from the file
  for (unsigned i=0;i<_n1;i++) {        // logx
    for (unsigned j=0;j<_n2;j++) {      // logkt2

        double logx, logkt2, xg;
        indata >> logx >> logkt2 >> xg;
        // saving values in vectors could be improved here
        vec_logx[i] = logx;
        vec_logkt2[j] = logkt2;
        *(matrix_xg+i*_n2+j) = xg;
    }
  }
  indata.close();

  _norm = &_norm_F;
  _updf = new Interpolation2D(vec_logx, vec_logkt2,  matrix_xg, _int_type);
}

//------------------------------------------------------------------------------
// set gluon from 3d grid as updf
void UnintegratedGluon::_set_from_3dgrid(string filename) {

  // read the grid input file
  ifstream indata;
  indata.open(filename.c_str());
  if(!indata) { cerr << "Error: updf grid file could not be opened" << endl;
     exit(1);}

  // vectors for passing argument values
  vector<double> vec_logx(_n1);
  vector<double> vec_logkt2(_n2);
  vector<double> vec_logmu2(_n3);
  vector<double> vec_xg(_n3);


  // 3d array for passing function values at grid points
  double *matrix_xg = new double[_n1*_n2*_n3];


  // Read the grid from the file
  for (unsigned i=0;i<_n1;i++) {        // logx
    for (unsigned j=0;j<_n2;j++) {      // logkt2
      for (unsigned k=0;k<_n3;k++) {    // logmu2
        double logx, logkt2, logmu2, xg;
        indata >> logx >> logkt2 >> logmu2 >> xg;
        // saving values in vectors could be improved here
        vec_logx[i] = logx;
        vec_logkt2[j] = logkt2;
        vec_logmu2[k] = logmu2;
        *(matrix_xg+i*(_n2*_n3)+j*_n3+k) = xg;
	//cout << logx << " " << logkt2 << " " << logmu2 << " " << xg << endl;
      }
    }
  }
  indata.close();

  _norm = &_norm_F;
  _updf = new CubicSpline3D(vec_logx, vec_logkt2,  vec_logmu2, matrix_xg);
}

// KMS and BK-like normalization
double UnintegratedGluon::_norm_F(const vector<double> &x) {
  return 1.0;
}


//------------------------------------------------------------------------------
double UnintegratedGluon::xg(double logx, double logkt2, double logmu2) {
  vector<double> x(3);
  x[0] = logx;
  x[1] = logkt2;
  x[2] = logmu2;

  // take the value of interpolated gluon
  double val = _norm(x)*(_updf->interp(x));
  //cout << val << endl;

  // check if the corresponding x vector is in the domain of the grid
  // used for interpolation, if not set the gluon to zero
  if (_updf->eval_status() > 0) {
    *_errorfile << "Err001: xg(" << x[0] << ", " << x[1];
    if (_name == grid3d) *_errorfile << ", " << x[2];
    *_errorfile << ")" << endl;
    val=0.0;
  }
  
  // check if the gluon is positive, if not, set it to zero
  if (val < 0 ) {
    *_errorfile << "Err002: xg(" << x[0] << ", " << x[1];
    if (_name == grid3d) *_errorfile << ", " << x[2];
    *_errorfile << ") = " << val << endl;
    val=0.0;
  }

  return val;
}

//------------------------------------------------------------------------------
/*
UnintegratedGluonName UnintegratedGluon::_determine_grid_type(string filename) {

  if (filename.find("KSlinear")  != string::npos) {
    return grid2d;
  } else if (filename.find("KSnonlinear")  != string::npos) {
    return grid2d;
  } else if (filename.find("KShardscalelinear")  != string::npos) {
    return grid3d;
  } else if (filename.find("KShardscalenonlinear")  != string::npos) {
    return grid3d;
  } else if (filename.find("KSsudFqgnonlinear.dat")  != string::npos) {
    return grid3d;
  } else if (filename.find("KSsudFggnonlinear.dat")  != string::npos) {
    return grid3d;
  } else {
    cerr << "Error: Unknown grid file" << endl; exit(1);
  }

}
*/

