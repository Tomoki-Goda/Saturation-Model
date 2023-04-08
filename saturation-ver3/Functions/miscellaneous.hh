#ifndef MISCLELLANEOUS_HH
#define MISCLELLANEOUS_HH


inline double  modx(const double  x, const double  Q2, const  double  mf2){
#if MODX==1
	return( (x*(1+4*mf2/Q2)));
#elif MODX==0 
	return( x);
#endif
}
inline double alpha(const double mu2){
	return 4.0/(9.0 *log( ((mu2>2*LQCD2)?(mu2):(2.0*LQCD2))/LQCD2));
}

static double change_var(double & var,double &  jac,const double min, const double max,const double c){//This version is (in theory) regular at max->Inf
	double den=( (c==1)?(1):(c+var*(1-c)) );
	jac= ( (min==0.0)?(c*pow(den,-2)*max):(c*pow(den,-2)*(max-min)) ) ;
	var= ( (min==0.0)?(max*var):((max*var+c*min*(1-var)) ))/den;
	//var= (max*var+min*c*(1-var))/den;
	
//#if TEST==1	
	if(var>max) {
		if(fabs((var-max)/max)>1.0e-15){
			printf("value below limit %.3e -> %.3e [%.3e, %.3e] diff %.3e, c=%.3e\n",(1-den)/(1-c),var,min,max,var-max, c);
		}
		var=max;
	}else if(var<min){
		if(fabs((min-var)/min)>1.0e-15){
		printf("value below limit %.3e -> %.3e [%.3e, %.3e] diff %.3e, c=%.3e\n",(1-den)/(1-c),var,min,max,min-var, c);
		}
		var=min;
	}
//#endif
	return var;
}
template<typename T>static double deriv(T & func,double y, double x,double h,int i) {
	//////////////////////////////////////////////////////////////////////////////////////////////
	//take ith numerical derivative with step h wrt the second argument "x" while the first arg 'y' is fixed 
	//////////////////////////////////////////////////////////////////////////////////////////
	double xi[5];
	double c2[3];//={-1.0/12,4.0/3,-5.0/2,4.0/3,-1.0/12};
	int c;
	switch(i){
	case 1:
		c=-1;
		c2[0]= 1.0/12.0;
		c2[1]=-2.0/3.0;
		c2[2]= 0;
		break;
	case 2:
		c=1;
		c2[0]=-1.0/12.0;
		c2[1]= 4.0/3.0;
		c2[2]=-5.0/2.0;
	     	break;
	default:
		printf("deriv:: Unknown option %d\n",i);
		getchar();
	}

	for(int j=0;j<5;j++){
		xi[j]=x+(j-2)*h;
	}
	double val=0;
	val =c2[0]*( func(y,xi[0])+c*func(y,xi[4]) );
	val+=c2[1]*( func(y,xi[1])+c*func(y,xi[3]) );
	if(i!=1){
		val+=c2[2]*  func(y,xi[2]);
	}
	val*=pow(h,-i);
	return(val);	
}
template<typename T>static double deriv2(T & func,double y, double x,double h,int i) {
	//////////////////////////////////////////////////////////////////////////////////////////////
	//take ith numerical derivative with step h wrt the second argument "x" while the first arg 'y' is fixed 
	//////////////////////////////////////////////////////////////////////////////////////////
	double xi[3];
	double c2[2];//={-1.0/12,4.0/3,-5.0/2,4.0/3,-1.0/12};
	int c;
	switch(i){
	case 1:
		c=-1;
		c2[0]= -1.0/2.0;
		c2[1]= 0;
		break;
	case 2:
		c=1;
		c2[0]= 1.0;
		c2[1]=-2.0;
	     	break;
	default:
		printf("deriv:: Unknown option %d\n",i);
		getchar();
	}

	for(int j=0;j<3;j++){
		xi[j]=x+(j-2)*h;
	}
	double val=0;
	val =c2[0]*( func(y,xi[0])+c*func(y,xi[2]) );
	if(i!=1){
		val+=c2[1]*  func(y,xi[1]);
	}
	val*=pow(h,-i);
	return(val);	
}
inline double min(double a,double b){
	return((a>b)?(b):(a));
}
		
inline double max(double a,double b){
	return((a>b)?(b):(a));
}
#endif
