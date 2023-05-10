	#include<cmath>
#include<iostream>
#include<vector>
#include<string>
#include<complex>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

/////////////////////////////////////////////
//arguments are  
//		directory/Name of gluon grid file
//		directory/Name of sigma grid file
//		directory to save the result
////////////////////////////////////////////

int main(int c , char** argv){
		//std::string dir(argv[2]);
		char filename[500];
		
		unsigned k2_npts=500, x_npts=500;
		double *k2_array,*x_array,*sigma_array;
		gsl_interp_accel *  k2_accel_ptr,*x_accel_ptr;
		gsl_spline2d *  spline_ptr;
	
		x_array=(double*)malloc(x_npts*sizeof(double));
		k2_array=(double*)malloc(k2_npts*sizeof(double));
		sigma_array=(double*)malloc(x_npts*k2_npts*sizeof(double));
		
	
		x_accel_ptr = gsl_interp_accel_alloc ();
		k2_accel_ptr = gsl_interp_accel_alloc ();
		
		FILE* file=fopen(argv[1],"r");
		
		double x,k2,val,xprev=0.0;
		double x0;
		for(int j=0;j<x_npts;++j){
			for(int i=0;i<k2_npts;++i){
				fscanf(file,"%lf %lf %lf",&x,&k2,&val);
				//printf("%lf %lf %lf %d\n",x,k2,val,feof(file));
				if(i==0){
					x_array[j]=exp(x);
				}
				if(j==0&&i==0){
					xprev=x;
				}
				
				if(x==xprev){
					if(j==0){
						k2_array[i]=exp(k2);
						//printf("r: %.3e\n",k2_array[i]);
					}
				}else if(j==0){
					k2_npts=i;
				}
				sigma_array[j*k2_npts+i]=val;		
			}
		
			//printf("%d %d\n",x_npts,k2_npts);
			if(feof(file)){
				x_npts=j+1;
				break;
			}
			//printf("x: %.3e\n",x_array[j]);
		}
		printf("Gluon grid shape %d %d\n",x_npts,k2_npts);
		spline_ptr = gsl_spline2d_alloc(gsl_interp2d_bicubic,k2_npts, x_npts);
		gsl_spline2d_init (spline_ptr,k2_array, x_array, sigma_array, k2_npts, x_npts);
		fclose(file);
		printf("GLUON\n");	
		int xplot[2]={2,6};
		for(int j=0;j<2;++j){
			sprintf(filename,"%s/%s%d.txt",argv[3],"gluon-",xplot[j]);
			file=fopen(filename,"w");
			x=pow(10,-xplot[j]);
			for(int i=0;i<100;++i){
				k2=1.0e-3*pow(1.0e+6,((double)i)/99);
				val=gsl_spline2d_eval(spline_ptr,k2,x,k2_accel_ptr, x_accel_ptr);
				fprintf(file ," %.5e %.5e\n", k2,val);
			}
		}	
		fclose(file);
///////////////////////////////////////////
//
		printf("SATURATION\n");
//
//
//////////////////////////////////////////
		double diff=1.0e-1;
		double valprev=0;
		int flag=0;
		int counter=0;
		double grad=0;
		sprintf(filename,"%s/%s.txt",argv[3],"saturation");
		file=fopen(filename,"w");
		for(int i=0;i<25;++i){
			k2=0.01;
			diff=1.0e-1;
			valprev=0;
			flag=0;
			counter=0;
			grad=0;
			x=2.5e-8*pow(1.0e+6,((double)i)/24);
			while(counter++<250){
				val=gsl_spline2d_eval_deriv_x(spline_ptr,k2, x,k2_accel_ptr, x_accel_ptr);
				if(fabs(val)<1.0e-5){
					break;
				}else if(val>0&&flag==0){
					k2+=diff;
					diff*=1.5;
				}else{
					flag=1;
					grad=(val-valprev)/diff;
					diff=-val/grad;
					k2+=diff;
				}
				valprev=val;
			}
			fprintf(file,"%.5e %.5e\n",x,k2);
			printf("%.5e %.5e\n",x,k2);
		}
		
		fclose(file);
/////////////////////////////////////////////////////////
//
		printf("SIGMA\n");
//
/////////////////////////////////////////////////////////		
		
		file=fopen(argv[2],"r");
		
		xprev=0.0;
		double r;
		for(int j=0;j<x_npts;++j){
			for(int i=0;i<k2_npts;++i){
				fscanf(file,"%lf %lf %lf",&x,&r,&val);
				//printf("%lf %lf %lf %d\n",x,k2,val,feof(file));
				if(i==0){
					x_array[j]=exp(x);
				}
				if(j==0&&i==0){
					xprev=x;
				}
				
				if(x==xprev){
					if(j==0){
						k2_array[i]=exp(r);
						printf("r: %.3e\n",k2_array[i]);
					}
				}else if(j==0){
					k2_npts=i;
				}
				sigma_array[j*k2_npts+i]=val;		
			}
		
			//printf("%d %d\n",x_npts,k2_npts);
			if(feof(file)){
				x_npts=j+1;
				break;
			}
			printf("x: %.3e\n",x_array[j]);
		}
		printf("SIGMA grid shape %d %d\n",x_npts,k2_npts);
		spline_ptr = gsl_spline2d_alloc(gsl_interp2d_bicubic,k2_npts, x_npts);
		gsl_spline2d_init (spline_ptr,k2_array, x_array, sigma_array, k2_npts, x_npts);
		fclose(file);
		
		//int xplot[2]={2,6};
		for(int j=0;j<2;++j){
			sprintf(filename,"%s/%s%d.txt",argv[3],"sigma-",xplot[j]);
			file=fopen(filename,"w");
			x=pow(10,-xplot[j]);
			for(int i=0;i<100;++i){
				r=1.0e-4*pow(5.0e+5,((double)i)/99);
				val=gsl_spline2d_eval(spline_ptr,r,x,k2_accel_ptr, x_accel_ptr);
				fprintf(file ," %.5e %.5e\n", r,val);
			}
		}	
		fclose(file);
		
		
		
		gsl_spline2d_free (spline_ptr);
		gsl_interp_accel_free (x_accel_ptr);
		gsl_interp_accel_free (k2_accel_ptr);
		free(k2_array);
		free(x_array);
		free(sigma_array);
		return 0;
		}	
