#include "InterpolationKS.hh"

//------------------------------------------------------------------------------
//
//   2D INTERPOLATION
//
//------------------------------------------------------------------------------

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
Interpolation2D::Interpolation2D(const vector<double> &x1vv, 
                             const vector<double> &x2vv, 
			     double* ym,
			     gsl_interp_type int_type)
                             :_int_type(int_type),
                               n1(x1vv.size()), n2(x2vv.size()), y(ym),
			       x1v(x1vv), x2v(x2vv) {
    
    acc.resize(n1);
    spline.resize(n1);


    // interpolate in x2 for all points of {x1}  set
    for (unsigned i=0;i<n1;i++) {  // x1
       acc[i] = gsl_interp_accel_alloc ();
       spline[i] = gsl_spline_alloc (&_int_type, n2);
       vector<double> y2;  
       for (unsigned j=0;j<n2;j++)
         y2.push_back(*(y + (i*n2) + j)); // x2

       gsl_spline_init (spline[i], &x2v[0], &y2[0], n2);
    }
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
Interpolation2D::~Interpolation2D() {
    // clean the memory corresponding to pointers from spline and acc matrices
    for (unsigned i=0;i<n1;i++) {
       gsl_interp_accel_free (acc[i]);
       gsl_spline_free (spline[i]);
    }
}

//-----------------------------------------------------------------------
vector<pair <double, double> > Interpolation2D::grid_limits() const {

  pair <double, double> x1_limits(*(x1v.begin()),*(x1v.rbegin()));
  pair <double, double> x2_limits(*(x2v.begin()),*(x2v.rbegin()));

  vector<pair <double, double> > limits;
  limits.push_back(x1_limits);
  limits.push_back(x2_limits);

  return limits;
}
  
//-----------------------------------------------------------------------
// return value from 2d spline interpolation
// get result for interpolation in x2 and use it to do interpolations in x1
//-----------------------------------------------------------------------
double Interpolation2D::interp (const vector<double> &x) {

    double x1=x[0], x2=x[1];
    double val;
    bool first_x2out = true;
    _eval_status = 0;
   
    // get values for {x1} grid
    double y1[n1];
    for (unsigned i=0;i<n1;i++) {  // x1
      int gsl_status = gsl_spline_eval_e (spline[i], x2, acc[i], &y1[i]);
      // if x2 is outside the domain, second bit is set to 1
      if (gsl_status == GSL_EDOM && !(_eval_status & 2)) { 
         _eval_status = _eval_status | 2; }
    }

    // interpolate in x1
    gsl_interp_accel* acc1d = gsl_interp_accel_alloc ();
    gsl_spline* spline1d = gsl_spline_alloc (&_int_type, n1);
    gsl_spline_init (spline1d, &x1v[0], y1, n1);

    // save the final value in 'val' and clean memory related to spline
    int gsl_status = gsl_spline_eval_e (spline1d, x1, acc1d, &val);
    // if x1 is outside the domain, first bit is set to 1
    if (gsl_status == GSL_EDOM) { _eval_status = _eval_status | 1; }
    gsl_interp_accel_free (acc1d);
    gsl_spline_free (spline1d);

    return  val;
}

//------------------------------------------------------------------------------
//
//   3D INTERPOLATION
//
//------------------------------------------------------------------------------

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
CubicSpline3D::CubicSpline3D(const vector<double> &x1vv, 
                             const vector<double> &x2vv, 
                             const vector<double> &x3vv, 
			     double* ym)
                             : n1(x1vv.size()), n2(x2vv.size()),n3(x3vv.size()),
			      y(ym), x1v(x1vv), x2v(x2vv), x3v(x3vv) {
    
    acc.resize(n1);
    spline.resize(n1);

    // interpolate in z for all points of (x1, x2) grid
    for (unsigned i=0;i<n1;i++) {    // x1
      for (unsigned j=0;j<n2;j++) {  // x2
        acc[i].push_back(gsl_interp_accel_alloc ());
        spline[i].push_back(gsl_spline_alloc (gsl_interp_linear, n3));
        vector<double> y3;  
        for (unsigned k=0;k<n3;k++) 
          y3.push_back(*(y +i*(n2*n3)+j*n3+k)); // x2
        gsl_spline_init (spline[i][j], &x3v[0], &y3[0], n3);
     }
    }

}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
CubicSpline3D::~CubicSpline3D() {
    // clean the memory corresponding to pointers from spline and acc matrices
    for (unsigned i=0;i<n1;i++) {
      for (unsigned j=0;j<n2;j++) {
       gsl_interp_accel_free (acc[i][j]);
       gsl_spline_free (spline[i][j]);
      }
    }
}

  
//-----------------------------------------------------------------------
/// return value from 3d spline interpolation
/// get result for interpolation in x3 and use it to do two consecutive
/// interpolations in x2 and x1
//-----------------------------------------------------------------------
double CubicSpline3D::interp (const vector<double> &x) {

    double x1=x[0], x2=x[1], x3=x[2];
    double val;
    _eval_status = 0;

    // get values for (x1,x2) grid
    double y2[n1][n2]; 
    for (unsigned i=0;i<n1;i++)   // x1
      for (unsigned j=0;j<n2;j++) { // x2 
        int gsl_status = gsl_spline_eval_e (spline[i][j], x3, acc[i][j], &y2[i][j]);
        // if x3 is outside the domain, third bit is set to 1
        if (gsl_status == GSL_EDOM && !(_eval_status & 4)){
           _eval_status = _eval_status | 4; }//cout << "Third bit" << endl; }
        
      }

    // interpolate in x2
    vector<gsl_interp_accel*> acc2d (n2);
    vector<gsl_spline*> spline2d (n2);
    for (unsigned i=0;i<n1;i++) { // loop over x1
      acc2d[i] = gsl_interp_accel_alloc();
      spline2d[i] = gsl_spline_alloc(gsl_interp_linear, n2);
      gsl_spline_init(spline2d[i], &x2v[0], &y2[i][0], n2);
    }

    // get values for (x1) grid
    double y1[n1];
    for (unsigned i=0;i<n1;i++) {
      int gsl_status = gsl_spline_eval_e (spline2d[i], x2, acc2d[i], &y1[i]);
     // if x2 is outside the domain, second bit is set to 1
     if (gsl_status == GSL_EDOM && !(_eval_status & 2)){
        _eval_status = _eval_status | 2; } //cout << "Second bit" << endl;}
    }

    // interpolate in x1
    gsl_interp_accel* acc1d = gsl_interp_accel_alloc ();
    gsl_spline* spline1d = gsl_spline_alloc (gsl_interp_linear, n1);
    gsl_spline_init (spline1d, &x1v[0], y1, n1);

    // save the final value in 'val' and clean memory related to spline
    int gsl_status = gsl_spline_eval_e (spline1d, x1, acc1d, &val);
    // if x1 is outside the domain, first bit is set to 1
    if (gsl_status == GSL_EDOM) { _eval_status = _eval_status | 1; }//cout << "First bit" << endl;}   

    for (unsigned i=0;i<n1;i++) {
      gsl_interp_accel_free (acc2d[i]);
      gsl_spline_free (spline2d[i]);
    }
    gsl_interp_accel_free (acc1d);
    gsl_spline_free (spline1d);

    //cout << "x=" << x1 << " kt=" << x2 << " mu=" << x3 << " xg= " << val << endl;
    return  val;
}
//-----------------------------------------------------------------------
vector<pair <double, double> > CubicSpline3D::grid_limits() const {


  pair <double, double> x1_limits(*(x1v.begin()),*(x1v.rbegin()));
  pair <double, double> x2_limits(*(x2v.begin()),*(x2v.rbegin()));
  pair <double, double> x3_limits(*(x3v.begin()),*(x3v.rbegin()));

  vector<pair <double, double> > limits;
  limits.push_back(x1_limits);
  limits.push_back(x2_limits);
  limits.push_back(x3_limits);

  return limits;
}
