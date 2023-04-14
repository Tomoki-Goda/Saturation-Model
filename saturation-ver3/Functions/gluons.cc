#include"gluons.hh"
///////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////

std::complex<double> Collinear_Gluon::gammatilde(const std::complex<double>& n)const{
	std::complex<double> n1,n2,n3,l1,l2,t1,t2,t3,cx,value;
	double m1,m2,rl;
	int k=0;
	n1 = n+1.0;
	n2 = n+2.0;
	n3 = n+3.0;

	l1 = conj(n)*conj(n1);
	m1 = abs(n)*abs(n)*abs(n1)*abs(n1);
	t1 = (1.0/m1)*l1;

	l2 = conj(n2)*conj(n3);
	m2 = abs(n2)*abs(n2)*abs(n3)*abs(n3);
	t2 = (1.0/m2)* l2;

	t3 = digamma(n2);

	cx = (t1+t2)-t3;
	rl = 11.0/2.0-NF/3.0-6.0*GAMMA_E;

	value = 6.0*cx + rl;
	return value;
}
//double Collinear_Gluon::integrand(const double y,const std::vector<double> &par)const {
double Collinear_Gluon::operator()(const double y,const std::vector<double>&par)const {
	std::complex<double> n0,n1,n2,g1,g2,gt,ex,l;
	double val;
	double m;

	const double Yg=par[0], tg=par[1];
	const double lambda_g=par[2];
	gsl_sf_result resr,resi;
	std::complex<double>comp(0,1);
	n0 = n_0+y*comp;
	n1 =-lambda_g+n_0+comp*y;
	n2 =-lambda_g+beta+n_0+comp*y;

	gsl_sf_lngamma_complex_e(n1.real(),n1.imag(),&resr,&resi );
	g1=exp(resr.val+comp*resi.val);
	gsl_sf_lngamma_complex_e(n2.real(),n2.imag(),&resr,&resi );
	g2=exp(resr.val+comp*resi.val);

	gt = tg *gammatilde(n0 );
	ex = exp(comp*y* Yg+gt);
	l = g1*conj(g2);
	m = abs(g2)*abs(g2);
	val = ((1.0/m)*l*ex).real();
	if(not(std::isfinite(val))){
		return 0;
	}
	return val;
}
//////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////
double Collinear_Gluon::operator()(const double x, const double QQ,const double A_g,const double l_g)const  {
	
	double value;
	const double bprim = 33.0/6.0-NF/3.0;
	const std::vector<double> par={
		log(1/x),
		(1/bprim)*log(log(QQ/LQCD2)/log(Q0/LQCD2)),
		l_g
	};
	const double normalization=A_g*exp(n_0* par[0] )*dgammafbeta;
	value=dclenshaw<const Collinear_Gluon&,const std::vector<double>&>(cc,*this,par,0,150,1.0e-10,1.0e-15);  
	value=normalization*value;
	if(!std::isfinite(value)){
		return 0;
	}
	return value ;
}
////////////////////////////////////////////////////////////////////////
// 1D Approximation at fixed x
////////////////////////////////////////////////////////////////////////
inline double Chebyshev1D_Collinear_Gluon::operator()(double *arg,const Collinear_Gluon& xg)const{
	double q2=change_var_revert_log(q2min,q2max, arg[0]);//q2min*pow(q2max/q2min,arg[1]);
	return( xg(*fixx,q2,A_g,l_g) );	
}

int Chebyshev1D_Collinear_Gluon::init(double q2min, double q2max, double A_g,double l_g ){
	this->q2min=q2min;
	this->q2max=q2max;
	this->A_g=A_g;
	this->l_g=l_g;
	return 0;
}
void Chebyshev1D_Collinear_Gluon::set_x(const double &x){
	fixx=&x;
	//printf("Chebyshev1D_Collinear_Gluon x set to %.2e",x);
	cheb_coeff<Chebyshev1D_Collinear_Gluon,const Collinear_Gluon&>(cheb[0],*this,xg);
	double val1, val2, mu2;
	//printf("nf= %d\n", NF);
	
////////////////  TEST  ////////////
 /*  	for(int i=0;i<5;++i){
		mu2=2*q2min*pow(q2max/(4*q2min),((double)i)/4);
		val1=xg(x,mu2,A_g,l_g);
		val2=(*this)(x,mu2,A_g,l_g);
		if(fabs((val1-val2)/(val1+val2))>1.0e-5&&fabs(val1-val2)>1.0e-7){
			printf("Chebyshev1D_Collinear_Gluon:: val= %.2e %.2e, diff= %.2e, at x=%.2e Q2=%.2e, Ag=%.2e lg=%.2e\n",val1,val2,val1-val2,x,mu2,A_g, l_g );
		}
	}
*/	
///////////////////////////////////
}

double Chebyshev1D_Collinear_Gluon::operator()(const double x,const double Q2, double A_g,double l_g )const{
	if(x!=*fixx){
		printf("Chebyshev1D_Collinear_Gluon:: Error: x does not match. input=%.3e internal x= %.3e diff = %.3e\n",x,*fixx, x-*fixx);
	}
	if(Q2>q2max){
		printf("Q2 too large: %.3e < %.3e < %.3e\n",q2min,Q2,q2max );
	}
	if(Q2<q2min){
		printf("Q2 too small: %.3e < %.3e < %.3e\n",q2min,Q2,q2max );
	}
	if((this->A_g!=A_g)||(this->l_g!=l_g)){
		printf("Chebyshev1D_Collinear_Gluon:: Parameters don't match %.2e %.2e %.2e %.2e\n",this->A_g,A_g,this->l_g,l_g );
	}
	double res;
	double arg[]={
		change_var_compactify_log(q2min,q2max,Q2 )
	};
	
	res=chebyshev(cheb[0],arg);
	//res=A_g*pow(x,-l_g);
	return(res);
}


//////////////////////////////////////////////////////////////////////
// 2D approximation of xg(x,Q^2)
//////////////////////////////////////////////////////////////////////
inline int Chebyshev_Collinear_Gluon::init(double xmin,double xmax,double q2min, double q2max, double A_g,double l_g ){
	this->xmin=xmin;
	this->xmax=xmax;
	this->q2min=q2min;
	this->q2max=q2max;
	this->A_g=A_g;
	this->l_g=l_g;
	cheb_coeff<Chebyshev_Collinear_Gluon,const Collinear_Gluon&>(cheb[0],*this,xg);
	return 0;
}
inline void Chebyshev_Collinear_Gluon::set_x(const double &x){
	fixx=&x;
	if(x>0){
		chebyshev_reduce(cheb[0], cheb[1],change_var_compactify_log(xmin,xmax,x ), 0 );
	}
}
inline double Chebyshev_Collinear_Gluon::operator()(double *arg,const Collinear_Gluon& xg)const{
	double x=change_var_revert_log(xmin,xmax, arg[0]);//xmin*pow(xmax/xmin,arg[0]);
	double q2=change_var_revert_log(q2min,q2max, arg[1]);//q2min*pow(q2max/q2min,arg[1]);
	return( xg(x,q2,A_g,l_g) );	
}
double Chebyshev_Collinear_Gluon::operator()(const double x,const double Q2){

	if(x>xmax){
		printf("x too large: %.3e < %.3e < %.3e\n",xmin,x,xmax );
	}
	if(x<xmin){
		printf("x too small: %.3e < %.3e < %.3e\n",xmin,x,xmax );
	}
	if(Q2>q2max){
		printf("Q2 too large: %.3e < %.3e < %.3e\n",q2min,Q2,q2max );
	}
	if(Q2<q2min){
		printf("Q2 too small: %.3e < %.3e < %.3e\n",q2min,Q2,q2max );
	}
	double res;
	///////////////////////////////////////////////////////////////////////////////////////////		
	// This block should be removed if x varies often.
	// Only useful if x is fixed and Q changes rapidly.
	///////////////////////////////////////////////////////////////////////////////////////////
	if(fixx!=NULL&&*fixx>0){
		if(x!=*fixx){
			printf("Chebyshev1D_Collinear_Gluon:: Error: x does not match. input=%.3e internal x= %.3e diff = %.3e\n",x,*fixx, x-*fixx);
		}
		double arg[]={
				change_var_compactify_log(q2min,q2max,Q2 )
		};
		res=chebyshev(cheb[1],arg);

	}else{
		double arg[]={
			change_var_compactify_log(xmin,xmax,x ),
			change_var_compactify_log(q2min,q2max,Q2 )
		};
		res=chebyshev(cheb[0],arg);
	}	
	return(res);
}	



