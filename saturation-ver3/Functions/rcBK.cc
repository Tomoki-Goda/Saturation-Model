#include<cmath>
#include<iostream>
#include<vector>
#include<string>
#include<complex>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include"constants.h"
#include"Kahn.hh"
#include"clenshaw.hh"
#include <gsl/gsl_sf.h>
#include"Levin.hh"
#define SUDAKOV 0
#define SECTOR_MAX 100


class Sigma_rcBK_integrand{
		//double *r_array=NULL,*sigma_array=NULL;
		const unsigned r_npts=160, x_npts=101;
		double *r_array,*x_array,*sigma_array;
		gsl_interp_accel *  r_accel_ptr,*x_accel_ptr;
		gsl_spline2d *  spline_ptr;
	public:
		Sigma_rcBK_integrand(){
			x_array=(double*)malloc(x_npts*sizeof(double));
			r_array=(double*)malloc(r_npts*sizeof(double));
			sigma_array=(double*)malloc(x_npts*r_npts*sizeof(double));
		}
		~Sigma_rcBK_integrand(){
			gsl_spline2d_free (spline_ptr);
			gsl_interp_accel_free (x_accel_ptr);
			gsl_interp_accel_free (r_accel_ptr);
			free(r_array);
			free(x_array);
			free(sigma_array);
		}
		
		int init(std::string name){
			x_accel_ptr = gsl_interp_accel_alloc ();
			r_accel_ptr = gsl_interp_accel_alloc ();
			spline_ptr = gsl_spline2d_alloc(gsl_interp2d_bicubic,r_npts, x_npts);
			FILE* file=fopen(name.c_str(),"r");
			
			double x,r,val;
			double x0;
			for(int j=0;j<x_npts;++j){
				for(int i=0;i<r_npts;++i){
					fscanf(file,"%lf %lf %lf",&x,&r,&val);
					sigma_array[(x_npts-j-1)*r_npts+i]=val;
					if(j==0){
						r_array[i]=sqrt(r);
						printf("r: %.3e\n",r_array[i]);
					}
					
					if(i>0&&x0!=x){
						printf("check grid shape\n");
						exit(1);
					}else{
						x0=x;
					}					
				}
				x_array[(x_npts-j-1)]=x;
				//printf("x: %.3e\n",x);
			}
			gsl_spline2d_init (spline_ptr,r_array, x_array, sigma_array, r_npts, x_npts);
			return 0;
		}		
		double operator()(double r,const std::vector<double>& par){
			double kt=sqrt(par[0]),x=par[1];
			double val;		
			val=gsl_spline2d_eval(spline_ptr,r, x,r_accel_ptr, x_accel_ptr);
			val=(1-pow(1-val,9/4));
			val*=std::cyl_bessel_j(0,r *kt );
			val/=r;
			
			return val;
		}
		
};
class rcBK_Gluon{
		Sigma_rcBK_integrand integrand;
		CCIntegral cc=CCprepare(64,"dipole",1,5);
		double INT_PREC=1.0e-4;
	public:
		rcBK_Gluon(){
			integrand.init("./resBKlingrid.dat");
		}	
		
		double operator()(double x,double kt2,double mu2){
			Kahn accum=Kahn_init(3);
			const std::vector<double> par{kt2,x,mu2};
			double rmax=9.0e+1,rmin=2.0e-3;
			double sum_accel, err;
			double sum = 0;
			double val=0,val1=0,val2=0;
			Kahn_clear(accum);

			double scale=(PI)/sqrt(kt2);
			double imin=rmin;
			int sectors=(int)(rmax/scale);

			if(sectors>SECTOR_MAX||sectors<1||!std::isfinite(sectors)){
				sectors=SECTOR_MAX;
			}

			double imax=PI/(sqrt(kt2)*4); //forJ0 integral, this is efficient
						  
			int flag=0,pass=5,accel_len=6,accel_min=3;
			Levin lev(accel_len+2);
			while(imax+scale<imin){
				imax+=scale;
			}
			
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
				
				val=dclenshaw<Sigma_rcBK_integrand,const std::vector<double>&>(cc,integrand,par,imin,imax,pow(INT_PREC,2),10e-14);
				
				imin=imax;
				if(val==0.0||sum+val==sum){
					sectors=i;
					break;
					//i--;
					//continue;
				}
				lev.add_term(val);
				sum=lev.sum(i);
				
				
				if(i>=accel_min*accel_len&&(flag>1 || (3*(i/3))==i )){
					
					val1=lev.accel(i-accel_len,accel_len);
					if(flag>=1){
						if(fabs(2*(val1-val2)/(val1+val2))<INT_PREC/10||fabs(val2-val1)<pow(INT_PREC/5,2) ){
							//++flag;
							if(++flag==pass){
								sectors=i+1;
								break;
							}
						}else{
							flag=1;//reset
						}
					}else if(flag==0){
						flag=1;
					}
					
					val2=val1;
				}
			}
			--sectors;

			//if(sectors>=25){
			if(flag==pass){
				val1=lev.accel(sectors-accel_len,accel_len);
				val2=lev.accel(sectors-1-accel_len,accel_len);	
				//val=val1;
			}else if(sectors>=accel_min*accel_len){
				val1=lev.accel(sectors-accel_len,accel_len);
				val2=lev.accel(sectors-1-accel_len,accel_len);
				if(fabs(2*(val1-val2)/(val1+val2))>INT_PREC/2&&fabs(val2-val1)>2*pow(INT_PREC/5,2) ){
					printf("\n3: sum=%.3e \t lev=%.1e %.1e\t diff= %.2e\t %d rmax= %.1e x=%.1e kt2=%.1e last term=%.1e, INT_PREC=%.1e\n\n",
					sum,val1,val2,fabs(val1-val2),sectors,imax,x,kt2,val,INT_PREC);
				}
				//val=val1;
			}else{
				val1=sum;
			}
			val=0;
			val2=0;//just safety



			double diff=0;
			Kahn_free(accum);
			if(!std::isfinite(val)){
				printf("Dipole_Gluon:: % encountered 0 returned\n",val);
				val1=0;
			}
		//	#if WW==1
			val1*=2.0/(3.0*pow(PI,3));
			val1/=0.2;
			//val=((val<1.0e-5)?(1.0e-5):(val));
		//	#else
		//	val1*=3.0/(8.0*pow(PI,2));
		//	#endif			
			//printf("3:% e\n",val);
			return (val1);
		}

};

int main(){
	double X_MIN=5.0e-7,X_MAX=1.0e-2;
	double KT2_MIN=1.0e-2,KT2_MAX=1.0e+6;
	double &kt2max=KT2_MAX;
	rcBK_Gluon gluon;
	double val=0,x=0;
	double k2,mu2;
	double arr[7000]={0};
	FILE* outfile=fopen("rcBK-WW.dat","w");
	for(int k=0; k<100;++k){
		x=X_MIN*pow(X_MAX/X_MIN,((double)k)/99);

			
#pragma omp parallel private(k2,mu2)
{
			//double k2,mu2=0;
#pragma omp for schedule(dynamic)
			for(int j=0;j<100;++j){
				k2=KT2_MIN*pow(kt2max/KT2_MIN,((double)j)/99);
				printf("kt2= %.2e\n",k2);
/*#if SUDAKOV>=1
				for(int i=0;i<70;++i){
					mu2=1.0e-1*pow(1.0e+5/1.0e-1,((double)i)/69);
					arr[j*70+i]=gluon(x,k2,mu2);
				}
#else*/
				mu2=0;
				arr[j]=gluon(x,k2,mu2);
//#endif
				printf("\033[1A\033[2K\r");

			}
}
			for(int j=0;j<100;++j){
				k2=KT2_MIN*pow(kt2max/KT2_MIN,((double)j)/99);
/*#if SUDAKOV>=1
				for(int i=0;i<70;++i){
					mu2=1.0e-1*pow(1.0e+5/1.0e-1,((double)i)/69);
					val=arr[j*70+i];
					fprintf(outfile ,"%.10e\t%.10e\t%.10e\t%.10e\n",log(x),log(k2),log(mu2), val );
					//fprintf(outfile ,"%.10e\t%.10e\t%.10e\n",log(x),log(k2), val );
				}
#else*/	
				val=arr[j]*32.895/0.389 ;
				fprintf(outfile ,"%.10e\t%.10e\t%.10e\n",log(x),log(k2),val );
//#endif
				

			}
			//printf("\033[1A\033[2K\r");
			//gluon.set_max(kt2max,mu2);	
			
			//outfile=fopen(filenames,"a");
			//gluon.export_grid(outfile);
			//fclose(outfile);
		}
		fclose(outfile);
	
}

