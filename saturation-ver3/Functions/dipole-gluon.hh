

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

double change_var(double & var,double &  jac,const double min, const double max,const double c){//This version is (in theory) regular at max->Inf
	double den=( (c==1)?(1):(c+var*(1-c)) );
	jac= ( (min==0.0)?(c*pow(den,-2)*max):(c*pow(den,-2)*(max-min)) ) ;
	var= ( (min==0.0)?(max*var):((max*var+c*min*(1-var)) ))/den;
	//var= (max*var+min*c*(1-var))/den;
	
#if TEST==1	
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
#endif
	return var;
}

////////////////////////////////////////////////////////////////////
//
//GBW gluon
//
////////////////////////////////////////////////////////////////////
class Gluon_GBW{
	double sigma_0=0,lambda=0,x_0=0,mu02=0;
	//double x=0;
	std::string key;

	
	public:
		explicit Gluon_GBW(){
		}
		void init(const double *par){
				sigma_0 =(double)par[0];
				lambda	=(double)par[1];
				x_0	=(double)par[2];
#if MU02==0
				mu02 = par[3];
#else 
				mu02 = MU02;
#endif
		}
		~Gluon_GBW(){}
		
	public:

		//void set_x(double x){
		//	this->x=x;
		//}
		double operator()(const double x,const double k2,double mu2){
			if(x_0<1.0e-5||x_0>1.0e-3){
				return 0;
			}
			if(lambda<0.05||lambda>0.95){
				return 0;
			}
			double Qs2=pow(x_0/x,lambda);
			double val=3.0/(4*PI*PI)*k2/Qs2*exp(-k2/Qs2);
			if(std::isnan(val)==1){
				return(0);
			}
#if ALPHA_RUN==1
			val*=alpha(mu2+mu02)/0.2;
			//printf("%.2e %.2e\n",mu2,alpha(mu2));
#endif
#if THRESHOLD==1 
			double thresh_power=7;
			val*=pow(1-x,thresh_power);
#endif
			return (sigma_0*val) ;
		}
};


/////////////////////////////////////////////////////////////////////
//
//Dipole gluon
//
/////////////////////////////////////////////////////////////////////
//template <typename LSigma>class Dipole_Gluon{
class Dipole_Gluon{
		const double *par;
		Laplacian_Sigma integrand;
		CCIntegral cc=CCprepare(128,"dipole",4,2);
		double x;	

	public: 
		Dipole_Gluon(){
		}
		~Dipole_Gluon(){
			
		}
		inline void init(const int n,const double *par ){
			this->par=par;
#if LAPLACIAN==0
			integrand.init(n,par,'l');	
#elif LAPLACIAN==1
			integrand.init(n,par,'s');
#endif	
		}
		void set_x(double x){
			this->x=x;
			integrand.set_kinem(x);
			//integrand.approximate_thread(x);
		}
		double operator()(const double kt2,const double mu2){
			Kahn accum=Kahn_init(3);
			double rmax=R_MAX,rmin=R_MIN;
			const std::vector<double> par{kt2,x};
			double val=0;
			Kahn_clear(accum);
#if GBW_APPROX==1
			if(x>0.7){
				double qs2=(4*PI*PI*alpha(par[4])*par[1]*pow(x,-par[2])*pow(1-x,5.6))/(3*par[0]);
				val=2*par[0]*kt2/qs2*exp(-kt2/qs2);
				return(3.0/(8*PI*PI)*val);
			}
#endif
			const double minmax=50;
#if ADD_END>=0
			rmax=minmax;
#if MODEL==1
			rmax/=pow(1-x,2.5);
#endif
			if(rmax>R_MAX||!std::isfinite(rmax)){
				//printf("rmax %.3e reduced to %.3e\n",rmax,R_MAX );
				rmax=R_MAX;
			}
			const double scale=(2*PI)/sqrt(kt2);
			double imin=rmin;
			int j=(int)(rmax/scale);
			if(j<3){
				j=3;
			}
			
//#if IBP==0||IBP==2 //actually important. better cancellation for convergence.
			double imax=PI/(sqrt(kt2)*4);
			//double imax=PI/(sqrt(kt2)*2);
			//with this differences coming from the first term of ibp is small.
//#elif IBP==1
			//double imax=3*PI/(sqrt(kt2)*4);
//#endif
			
			int flag=0;
			//double val_prev=0;
			for(int i=0;i<j;++i){
				imax+=scale;
				if(imax>rmax){
					if(i>0){
						imax-=scale;
						j=i-1;
						break;
					}else{
						imax=rmax;
					}
				};
#if R_CHANGE_VAR==1
				val=dclenshaw<const Laplacian_Sigma,const std::vector<double>&>(cc,integrand,par,imin/(1+imin),imax/(1+imax),INT_PREC/10,INT_PREC/100);
#elif R_CHANGE_VAR==0
				val=dclenshaw<const Laplacian_Sigma,const std::vector<double>&>(cc,integrand,par,imin,imax,INT_PREC/10,INT_PREC/100);
#endif
				if(fabs(val)<1.0e-10){
					++flag;
					if(flag>6&&imax>minmax){
						break;//it is likely beyond this will be trivial
					}
				}else{
					flag=0;
				}
				accum+=val;
				imin=imax;
			}
			val=Kahn_total(accum);
#endif
			double diff=0;
//#if (IBP>=1&&ADD_END!=0)			
#if (IBP>=1)			
			diff+=integrand.constant(imax,par);
			diff-=integrand.constant(rmin,par);
			if(fabs(diff)>fabs(val/1.0e-9)&&fabs(diff)>1.0e-9){
				printf("inaccurat IBP val=%.1e diff=%.1e imax=%.2e rmax=%.2e scale= %.1e\n",val,diff,imax,rmax, scale);
				printf(" x= %.2e kt2=%.2e %.3e  %.3e\n",x,kt2, imax-(j*scale+3*PI/(sqrt(kt2)*4)) ,(imax*sqrt(kt2)-(PI/4))/PI);
			}
#if ADD_END!=0
			val+=diff;
#endif
#endif
			Kahn_free(accum);

			
			if(!std::isfinite(val)){
				val=0;
			}
			return (3.0/(8*PI*PI)*val);
		}
		
		
};

