/////////////////////////////////////////////////////////////////////
// Gluon_GBW af
// af.init(par) //double* par sigma parameters 
// af(x,kt2,mu2) // 
//
// Dipole gluon
// Dipole_Gluon af
// af.init(par,integ) //Integrand integ; see below
// af(x,kt2,mu2)
//
//WW_Gluon
////////////////////////////////////////////////////////////////////

//#if SIGMA_APPROX<2
//	typedef Gluon_Integrand Integrand;				
//#else
//	typedef Laplacian_Sigma Integrand;				
//#endif
//#include"series_sum.hh"
#include <gsl/gsl_sf.h>

#include"Levin.hh"

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

////////////////////////////////////////////////////////////////////
//
//GBW gluon
//
////////////////////////////////////////////////////////////////////
class Gluon_GBW{
	const double *sigpar;
	double sigma_0=0,lambda=0,x_0=0,mu02=0,thresh_power=0;
	double x=0;
	//std::string key;

	
	public:
		//explicit Gluon_GBW(const Gluon_GBW& rhs){
		//	sigpar=rhs.sigpar;
		//	init(sigpar);
		//	x=rhs.x;			
		//}
		
		explicit Gluon_GBW(){
		}
		void init(const double *par){
		 int count=0;
		 sigpar=par;
				sigma_0 =(double)par[count++];
				lambda	=(double)par[count++];
				x_0	=(double)par[count++];
#if MU02==0
				mu02 = par[count++];
#else 

				mu02 = MU02;
#endif
#if THRESHOLD==-1
				thresh_power=par[count++];
#else
				thresh_power=THRESHOLD;
#endif
		}
		~Gluon_GBW(){}
		
	public:
		//inline double operator()(const double x,const double k2,double mu2){
		//	set_x(x);
		//	return((*this)(k2,mu2));
		//}
		//void set_x(double x){
		//	this->x=x;
		//}
		double operator()(const double  x,const double k2,double mu2){
			if(x_0<1.0e-5||x_0>1.0e-3){
				return 0;
			}
			if(lambda<0.05||lambda>0.95){
				return 0;
			}
			double Qs2=pow(x_0/x,lambda);
#if THRESHOLD==-2
			Qs2*=pow(1-x,5);
#endif		
#if WW==1
			Qs2*=9.0/4.0;
			gsl_sf_result result;
			gsl_sf_gamma_inc_e(1.0e-20, k2/Qs2,&result);
			double val=result.val;
			val/=(3*pow(PI,3));
			val/=0.2;
#else			
			double val=3.0/(4*PI*PI)*k2/Qs2*exp(-k2/Qs2);
#endif
			if(std::isnan(val)==1){
				return(0);
			}
#if ALPHA_RUN==1
			val*=alpha(mu2+mu02)/0.2;
			//printf("%.2e %.2e\n",mu2,alpha(mu2));
#endif
#if THRESHOLD>0||THRESHOLD==-1 
			//double thresh_power=THRESHOLD;
			val*=pow(1-x,thresh_power);
#endif
			return (sigma_0*val) ;
		}
};
/*
inline double min(double a,double b){
	return((a>b)?(b):(a));
}
		
inline double max(double a,double b){
	return((a>b)?(b):(a));
}
*/
/////////////////////////////////////////////////////////////////////
//
//Dipole gluon
//
/////////////////////////////////////////////////////////////////////
/*//template <typename LSigma>class Dipole_Gluon{
//double eps(int n, int m, const double (&list)[SECTOR_MAX]){ 
double eps(int n, int m, const double *list){ 
//epsilon series convergence algorithm??
// I don't know the name but Seki-Aitken's generalization.
	//if(2*(m/2)!=m){
	//	printf("eps use multiple of 2, m=%d\n",m);
	//	m=2*((m+1)/2);
	//}
	if(m == -1){
		return(0);
	};
	if(m == 0){
	  	return(list[n]);
	};
	return(eps(n + 1, m - 2,list) + 1/(eps(n + 1, m - 1,list) - eps(n, m - 1,list)));
}*/
template<typename INTEG>class Dipole_Gluon{
		const double *par;
		//Laplacian_Sigma integrand;
		INTEG *integrand;
		CCIntegral cc=CCprepare(64,"dipole",1,5);
		
		
		//double fixx;

	public: 
		//Dipole_Gluon(const Dipole_Gluon&rhs ){
		//	par=rhs.par;
		//	integrand=rhs.integrand;
		//	x=rhs.x;
		//}
		Dipole_Gluon(INTEG& integ){
			integrand =& integ;
		}
		~Dipole_Gluon(){
			
		}
		inline void init(const double * const &par ){
			
			this->par=par;
			
			//init(n,par) if Laplacian_Sigma is used
//#if LAPLACIAN==0
//			integrand.init(par,'l');	
//#elif LAPLACIAN==1
//			integrand.init(par,'s');
//#endif	
		}
		void set_x(const double &x){
			//this->x=x;
			integrand->set_x(x);
		}
		double operator()(const double x,const double kt2,const double mu2){
			//Series_Sum ss(2);
			double sum_accel, err;
			double arr[SECTOR_MAX];
			double sum = 0;
			int n;
			//gsl_sum_levin_u_workspace * w = gsl_sum_levin_u_alloc (SECTOR_MAX);
			Kahn accum=Kahn_init(3);
			const std::vector<double> par{kt2,x,mu2};
			double rmax=R_MAX,rmin=R_MIN;
			double val=0,val1=0,val2=0;
			Kahn_clear(accum);
			
			double scale=(PI)/sqrt(kt2);
			double imin=rmin;
			int sectors=(int)(rmax/scale);
			
			if(sectors>SECTOR_MAX||sectors<1||!std::isfinite(sectors)){
				sectors=SECTOR_MAX;
			}
			
			/*if(sectors>10){
				scale*=2;
			}else{
				scale/=8;
				sectors*=8;
			}*/
			/*if(sectors>10){
				scale*=2;
			}*/
			/*
			if(sectors>20){
				scale*=2;
			}*/
#if IBP==1&&WW!=1
			double imax=3*PI/(sqrt(kt2)*4); //forJ0 integral, this is efficient
#else
			double imax=PI/(sqrt(kt2)*4); //forJ0 integral, this is efficient
#endif						      
						  
			int flag=0,pass=4,accel_len=6;
			Levin lev(SECTOR_MAX);
			for(int i=0;i<sectors;++i){
				imax+=scale;
				if(imax>rmax){
					if(i>0){
						imax-=scale;
						sectors=i;
						break;
					}else{
						imax=rmax;
					}
				};
				
#if R_CHANGE_VAR==1
				val=dclenshaw<INTEG,const std::vector<double>&>(cc,*integrand,par,imin/(1+imin),imax/(1+imax),pow(INT_PREC,2),10e-14);
#elif R_CHANGE_VAR==0
				val=dclenshaw<INTEG,const std::vector<double>&>(cc,*integrand,par,imin,imax,pow(INT_PREC,2),10e-14);
				//val=dgauss<INTEG,const std::vector<double>&>(*integrand,par,imin,imax,pow(INT_PREC,2),10e-14);
#endif
				//printf("%.3e\n",val);
				//getchar();
				if(val==0.0){
					sectors=i;
					break;
				}
				imin=imax;
				
				lev.add_term(val);
				sum=lev.sum(i);
				
				if(i>=25&&(flag>1 || (3*(i/3))==i )){
					//flag=0 untested
					//flag=1 tested without pass
					//flag>1 passed flag-1 times consecutively
					val1=lev.accel(i-accel_len,accel_len);
					if(flag>=1){
						if(fabs((val1-val2)/(val1+val2))<INT_PREC/5||fabs(val2-val1)<pow(INT_PREC/5,2) ){
							++flag;
						}else{
							if(flag>2){
								//printf("reset\n");
								//printf("0: sum=%.3e \t lev=%.1e %.1e %.3e\t %d rmax= %.2e x=%.2e kt2=%.2e last term=%.2e\n",sum,val1,val2,val1-val2,i,imax,x,kt2,val);
							}
							flag=1;//reset
						}

					}else if(flag==0){
						flag=1;
					}
					if(flag==pass){
						sectors=i+1;
						break;
					}
					val2=val1;
				}
			}
			--sectors;
			
			//if(sectors>=25){
			if(flag==pass){
				val1=lev.accel(sectors-accel_len,accel_len);
				val2=lev.accel(sectors-1-accel_len,accel_len);
				val=val1;
			}else if(sectors>=SECTOR_MAX/3&&sectors>accel_len){
				val1=lev.accel(sectors-accel_len,accel_len);
				val2=lev.accel(sectors-1-accel_len,accel_len);
				if(fabs(1-val2/val1)>INT_PREC&&fabs(val2-val1)>pow(INT_PREC,2) ){
					printf("2: sum=%.3e \t lev=%.1e %.1e\t diff= %.2e\t %d rmax= %.1e x=%.1e kt2=%.1e last term=%.1e\n",sum,val1,val2,fabs(val1-val2),sectors,imax,x,kt2,val);
				}
			}else{
				val=sum;
			}
			
			

			double diff=0;
#if (IBP>=1 && ADD_END!=0 && WW==0 )			
			diff+=integrand->constant(imax,par);
			if(fabs(diff)>fabs(val/1.0e-9)&&fabs(diff)>1.0e-9){
				printf("Dipole_Gluon:: inaccurat IBP val=%.1e diff=%.1e imax=%.2e rmax=%.2e scale= %.1e\n",val,diff,imax,rmax, scale);
				printf("Dipole_Gluon::  x= %.2e kt2=%.2e %.3e  %.3e\n",x,kt2, imax-(sectors*scale+3*PI/(sqrt(kt2)*4)) ,(imax*sqrt(kt2)-(PI/4))/PI);
			}
			val+=diff;
#endif//IBP
			Kahn_free(accum);
			if(!std::isfinite(val)){
				printf("Dipole_Gluon:: % encountered 0 returned\n",val);
				val=0;
			}
#if WW==1
			val*=2.0/(3.0*pow(PI,3));
			val/=0.2;
			//val=((val<1.0e-5)?(1.0e-5):(val));
#else
			val*=3.0/(8.0*pow(PI,2));
#endif			
			//printf("3:% e\n",val);
			return (val);
		}
		
		
};


 
/////////////////////////////////////////////////////////////////////
//
//WW gluon
//
/////////////////////////////////////////////////////////////////////
/*class WW_Gluon{
		const double *par;
		Integrand integrand;
		CCIntegral cc=CCprepare(64,"dipole",1,4);
		double x;	

	public: 
		WW_Gluon(const WW_Gluon&rhs ){
			par=rhs.par;
			integrand=rhs.integrand;
			x=rhs.x;
		}
		WW_Gluon(){
		}
		~WW_Gluon(){
			
		}
		inline void init(const int n,const double * const &par ){
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
			const std::vector<double> par{kt2,x};
			double rmax=R_MAX,rmin=R_MIN;
			double val=0;
			Kahn_clear(accum);

			const double minmax=50;
#if ADD_END>=0 //negative value is for testing purpose.
			rmax=minmax;
#if MODEL==1
			rmax/=pow(1-x,4);
#endif
			if(rmax>R_MAX||!std::isfinite(rmax)){
				rmax=R_MAX;
			}
			const double scale=(2*PI)/sqrt(kt2);
			double imin=rmin;
			int sectors=(int)(rmax/scale);
			if(sectors<3){
				sectors=3;
			}
			
			double imax=PI/(sqrt(kt2)*4); //forJ0 integral, this is efficient 
			int flag=0;
			//double val_prev=0;
			for(int i=0;i<sectors;++i){
				imax+=scale;
				if(imax>rmax){
					if(i>0){
						imax-=scale;
						sectors=i-1;
						break;
					}else{
						imax=rmax;
					}
				};
#if R_CHANGE_VAR==1
				val=dclenshaw<const Integrand,const std::vector<double>&>(cc,integrand,par,imin/(1+imin),imax/(1+imax),INT_PREC/(10*sectors),pow(INT_PREC,2));
#elif R_CHANGE_VAR==0
				val=dclenshaw<const Integrand,const std::vector<double>&>(cc,integrand,par,imin,imax,INT_PREC/(10*sectors),pow(INT_PREC,2));
#endif
				if(fabs(val+integrand.constant(imax,par))< pow(INT_PREC,2) ){
					++flag;
					if(flag>5&&imax>minmax){//if consecutively small 5 times
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
			//diff-=integrand.constant(rmin,par);
			if(fabs(diff)>fabs(val/1.0e-9)&&fabs(diff)>1.0e-9){
				printf("inaccurat IBP val=%.1e diff=%.1e imax=%.2e rmax=%.2e scale= %.1e\n",val,diff,imax,rmax, scale);
				printf(" x= %.2e kt2=%.2e %.3e  %.3e\n",x,kt2, imax-(sectors*scale+3*PI/(sqrt(kt2)*4)) ,(imax*sqrt(kt2)-(PI/4))/PI);
			}
#if ADD_END!=0
			val+=diff;
#endif//ADD_END
#endif//IBP
			Kahn_free(accum);
//Threshold 1-x^7 is in dipole sigma. 
			
			if(!std::isfinite(val)){
				val=0;
			}
			return (3.0/(8*PI*PI)*val);
		}
		
		
};


*/
