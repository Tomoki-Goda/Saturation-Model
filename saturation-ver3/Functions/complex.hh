/* 
 *    The complex numbers operations.
 *    File taken from Numerical Recipes.
 *    modified for c++
 * 
 */

typedef struct{
		double r=0;
		double i=0;
}  doublecomplex;

doublecomplex operator + (const doublecomplex & a,const doublecomplex & b ){
	doublecomplex c;
	c.i=a.i+b.i;
	c.r=a.r+b.r;
	return c;			
}
doublecomplex operator + (const double & a,const doublecomplex & b ){
	doublecomplex c;
	c.i=b.i;
	c.r=a+b.r;
	return c;			
}
doublecomplex operator + (const doublecomplex & b,const double & a ){
	doublecomplex c;
	c.i=b.i;
	c.r=a+b.r;
	return c;			
}

doublecomplex operator - (const doublecomplex & a,const doublecomplex & b ){
	doublecomplex c;
	c.i=a.i-b.i;
	c.r=a.r-b.r;
	return c;			
}
doublecomplex operator - (const double & a,const doublecomplex & b ){
	doublecomplex c;
	c.i=-b.i;
	c.r=a-b.r;
	return c;			
}
doublecomplex operator - (const doublecomplex & b,const double & a ){
	doublecomplex c;
	c.i=b.i;
	c.r=b.r-a;
	return c;			
}


doublecomplex operator * (const doublecomplex & a,const doublecomplex & b ){
	doublecomplex c;
	c.r=a.r *b.r -a.i*b.i;
	c.i=a.i *b.r +a.r*b.i;			
	return c;
}
doublecomplex operator * (const double & a,const doublecomplex & b ){
	doublecomplex c;
	c.r= a*b.r;
	c.i= a*b.i;			
	return c;
}
doublecomplex operator * (const doublecomplex & b,const double & a ) {
	doublecomplex c;
	c.r= a*b.r;
	c.i= a*b.i;			
	return c;
}



doublecomplex Complex(double re, double im){
	doublecomplex c;
	c.r=re;
	c.i=im;
	return c;
}
doublecomplex Conjg(const doublecomplex& z){
	doublecomplex c;
	c.r=z.r;
	c.i = -z.i;
	return c;
}
doublecomplex operator / (const doublecomplex & a,const double& b ){
	doublecomplex c;
	c.r=a.r/b;
	c.i=a.i/b;
			
	return c;
}
doublecomplex operator / (const doublecomplex & a,const doublecomplex & b ){
	doublecomplex c;
	double den=b.r*b.r+b.i*b.i;
	c=a*Conjg(b);
	c=c/den;
	return c;
}
double Cabs(const doublecomplex & b ){
	return sqrt(b.r*b.r+b.i*b.i);
}

doublecomplex Cexp(const doublecomplex& z) {

    doublecomplex c;
    c.r = exp(z.r)*cos(z.i);
    c.i = exp(z.r)*sin(z.i);
    return c;
}
doublecomplex Csqrt(const doublecomplex& z)
{
	doublecomplex c;
	double x,y,w,r;
	if ((z.r == 0.0) && (z.i == 0.0)) {
		c.r=0.0;
		c.i=0.0;
		return c;
	} else {
		x=fabs(z.r);
		y=fabs(z.i);
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.r >= 0.0) {
			c.r=w;
			c.i=z.i/(2.0*w);
		} else {
			c.i=(z.i >= 0) ? w : -w;
			c.r=z.i/(2.0*c.i);
		}
		return c;
	}
}

