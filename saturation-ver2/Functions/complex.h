/* 
 *    The complex numbers operations.
 *    File taken from Numerical Recipes.
 * 
 */


#include <math.h>

typedef struct {double r,i;} doublecomplex;


doublecomplex Cadd(doublecomplex a, doublecomplex b)
{
	doublecomplex c;
	c.r=a.r+b.r;
	c.i=a.i+b.i;
	return c;
}

doublecomplex Csub(doublecomplex a, doublecomplex b)
{
	doublecomplex c;
	c.r=a.r-b.r;
	c.i=a.i-b.i;
	return c;
}


doublecomplex Cmul(doublecomplex a, doublecomplex b)
{
	doublecomplex c;
	c.r=a.r*b.r-a.i*b.i;
	c.i=a.i*b.r+a.r*b.i;
	return c;
}

doublecomplex Complex(double re, double im)
{
	doublecomplex c;
	c.r=re;
	c.i=im;
	return c;
}

doublecomplex Conjg(doublecomplex z)
{
	doublecomplex c;
	c.r=z.r;
	c.i = -z.i;
	return c;
}

doublecomplex Cdiv(doublecomplex a, doublecomplex b)
{
	doublecomplex c;
	double r,den;
	if (fabs(b.r) >= fabs(b.i)) {
		r=b.i/b.r;
		den=b.r+r*b.i;
		c.r=(a.r+r*a.i)/den;
		c.i=(a.i-r*a.r)/den;
	} else {
		r=b.r/b.i;
		den=b.i+r*b.r;
		c.r=(a.r*r+a.i)/den;
		c.i=(a.i*r-a.r)/den;
	}
	return c;
}

double Cabs(doublecomplex z)
{
	double x,y,ans,temp;
	x=fabs(z.r);
	y=fabs(z.i);
	if (x == 0.0)
		ans=y;
	else if (y == 0.0)
		ans=x;
	else if (x > y) {
		temp=y/x;
		ans=x*sqrt(1.0+temp*temp);
	} else {
		temp=x/y;
		ans=y*sqrt(1.0+temp*temp);
	}
	return ans;
}

doublecomplex Csqrt(doublecomplex z)
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

doublecomplex RCmul(double x, doublecomplex a)
{
	doublecomplex c;
	c.r=x*a.r;
	c.i=x*a.i;
	return c;
}

/* complex exp*/
doublecomplex Cexp(doublecomplex z) {

    doublecomplex c;
    c.r = exp(z.r)*cos(z.i);
    c.i = exp(z.r)*sin(z.i);
    return c;
}
