#include"dipole-gluon.hh"
extern int INT_PREC;
void Gluon_GBW::init(const double *par){
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
	double Gluon_GBW::operator()(const double  x,const double k2,double mu2){
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
//////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////

void Dipole_Gluon::init(const double * const &par ){
	this->par=par;
}
void Dipole_Gluon::set_x(const double &x){
	integrand->set_x(x);
}

double Dipole_Gluon::operator()(const double x,const double kt2,const double mu2){
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

	#if IBP==1&&WW!=1
	double imax=3*PI/(sqrt(kt2)*4); //forJ0 integral, this is efficient
	#else
	double imax=PI/(sqrt(kt2)*4); //forJ0 integral, this is efficient
	#endif						      
				  
	int flag=0,pass=5,accel_len=8;
	Levin lev(10);
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
		val=dclenshaw<INTEG ,const std::vector<double>&>(cc,*integrand,par,imin/(1+imin),imax/(1+imax),pow(INT_PREC,2),10e-14);
	#elif R_CHANGE_VAR==0
		val=dclenshaw<INTEG ,const std::vector<double>&>(cc,*integrand,par,imin,imax,pow(INT_PREC,2),10e-14);
	#endif
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
		