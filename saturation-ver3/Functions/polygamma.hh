#include<iostream>
#include<complex>
#include<vector>
#include"control-default.h"
#include"constants.h"
typedef std::complex<double> dc;
static dc digamma(const dc& z){
	const double c[6]=
		{8.33333333333333333e-2, 
		-8.33333333333333333e-3,
		3.96825396825396825e-3,
		-4.16666666666666667e-3,
		7.57575757575757576e-3,
		-2.10927960927960928e-2};
	const double del=1.0e-15;
	double y=imag(z), x=real(z);
	dc v,u,h,r,p;
	dc comp(0,1);
	//u=z;
	if(fabs(x)<del&&fabs(x+lrint(x))<del){
		printf("error 2\n");
		exit(1);
	}
	if(x<0){
		 u=-z;
	}else{
		 u=z;
	}
	v=u;
	h=0;
	if(fabs(x)<15){
		h=1.0/v;
		for(int i=1;i<(15-int(fabs(x)));++i){
			v+=1.0;
			h+=1.0/v;
		}
		v+=1;
	}
	r=1.0/pow(v,2);
	p=c[5]*r;
	for(int i=4;i>=0;--i){
		p=r*(c[i]+p);
	}
	//h=sign[k]*(fct[k+1]*h+(v*(fct[k]+p)+0.5*fct[k+1])/pow(v,k1) );
	h=-(h+(v*p+0.5)/v );
	h+=log(v);
	if(x<0){
		v=PI*u;
		x=real(v);
		y=imag(v);
		double a=sin(x);
		double b=cos(x);
		double t=tanh(y);
		p=( b - (a*t*comp))/(a+( b*t*comp) );
		h+=1.0/u+PI*p;
	}
	return h;
}
/*
dc wpsipg(dc z,int k){
	const double c1 = pow(PI,2), c2 = 2*pow(PI,3), c3 = 2*pow(PI,4),
	      c4 = 8*pow(PI,5);
     	const double sign[]={-1,+1,-1,+1,-1};
       	const double fct[]= {0,1,1,2,6,24};
	const double c[5][6]={
		{8.33333333333333333e-2, 
		-8.33333333333333333e-3,
		3.96825396825396825e-3,
		-4.16666666666666667e-3,
		7.57575757575757576e-3,
		-2.10927960927960928e-2},
		{1.66666666666666667e-1,
		-3.33333333333333333e-2,
		2.38095238095238095e-2,
		-3.33333333333333333e-2,
		7.57575757575757576e-2,
		-2.53113553113553114e-1},
		{5.00000000000000000e-1,
		-1.66666666666666667e-1,
		1.66666666666666667e-1,
		-3.00000000000000000e-1,
		8.33333333333333333e-1,
		-3.29047619047619048e+0},
		{2.00000000000000000e+0,
		-1.00000000000000000e+0,
		1.33333333333333333e+0,
		-3.00000000000000000e+0,
		1.00000000000000000e+1,
		-4.60666666666666667e+1},
		{10,-7,12,-33,130,-691}};
	const double del=1.0e-15;

	double y=imag(z), x=real(z);
	dc u,v,h,r,p;
	dc comp(0,1);

	u=z;
	int k1;
	if (k<0||k>4){
		printf("error 1\n");
		exit(1);
	}	
	if(fabs(x)<del&&fabs(x+lrint(x))<del){
		printf("error 2\n");
		exit(1);
	}
	k1=k+1;
	
	if(x<0){
		u=-u;
	}
	v=u;
	h=0;
	if(fabs(x)<15){
		h=1.0/pow(v,k1);
		for(int i=1;i<(15-int(fabs(x)));++i){
			v+=1.0;
			h+=1.0/pow(v,k1);
		}
		v+=1;
	}
	r=1.0/pow(v,2);
	p=c[k][5]*r;
	for(int i=4;i>=0;--i){
		p=r*(c[k][i]+p);
	}
	h=sign[k]*(fct[k+1]*h+(v*(fct[k]+p)+0.5*fct[k+1])/pow(v,k1) );
	if(k==0){
		h+=log(v);
	}
	
	if(x<0){
		v=PI*u;
		x=real(v);
		y=imag(v);
		double a=sin(x);
		double b=cos(x);
		double t=tanh(y);
		p=( b - (a*t*comp))/(a+( b*t*comp) );
		if(k==0){
			h+=1.0/u+PI*p;
		}else if(k==1){
			h=-h+1.0/pow(u,2)+c1*(pow(p,2)+1.0);
		}else if(k==2){
			h=h+2.0/pow(u,3)+c2*(pow(p,2)+1.0);
		}else if(k==3){
			r=pow(p,2);
			h=-h+6.0/pow(u,4)+c3*((3.0*r+4.0)*r+1.0);
		}else if(k==4){
			r=pow(p,2);
			h=h+24.0/pow(u,5)+c4*p*((3.0*r+5.0)*r+2.0);
		}
	}
	return h;
}*/
