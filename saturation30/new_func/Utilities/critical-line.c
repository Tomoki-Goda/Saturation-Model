#include<stdio.h>
#include<string.h>
#include<math.h>
#include<stdlib.h>

//#include"control_tmp.h"
#include"control.h"
#include"../control-default.h"
#include"../constants.h"
 
#include"./plot.c"

extern double SIGMA(double ,double, double ,double*,double *);
extern void approx_xg(double*);
extern int parameter(double*,double*,double*);
		


int pick_interval(double(*func)(double , double **), double ** par , double *min ,double *max,double *differences, double aim){
	if (max<min){
		printf("max %f less than min%f \n", *max,*min);
		return 1;
	}
	int points_n=15;
	double val[points_n];
	double step = (*max-*min)/(points_n-1);
	double diff[points_n];

	double new_lim[2];
	//double diffs[2];
	double *diffs;
	//printf("min %.3e max %.e  step %.3e\n" ,*(min),*(max) ,step);
	for(int i= 0 ;i<points_n;i++){
		
		val[i]=(*func)(*min+i*step , par);
		diff[i]=(val[i]-aim);
		//if(i!=0 && (diff[i-1]*diff[i])/(pow(diff[i],2))<0){
		//printf("diff %.2e sigma %.2e aim %.2e \n", diff[i],val[i],aim);
		if(i!=0 && (diff[i-1]*diff[i])<0){
			//sign of difference changed -> went pass the point.
			new_lim[0]=*min+(i-1) *step;
			new_lim[1]=*min+i*step;
			diffs=(diff+i-1);
			break;
		}
		//printf("%d\t%d\n",(i-1),points_n);
		//printf("%d %.2e\t%.2e\n ",i, val[i],diff[i]);

		if((i+1)==points_n){
			printf("not found\n");
			for(int j=0;j<=i;j++){
				printf("at  %.2e\t sigma = %.2e\tdiff = %.2e\n ",*min+j*step ,val[j],diff[j]);
			}
			printf("\n");
			getchar();

		}
	}
	//printf("pick int %.3e\t%.3e\n",*diffs,*(diffs+1));

	
	if(new_lim[0]>new_lim[1]){
		*min=new_lim[1];
		*max=new_lim[0];
		differences[1]=diffs[1];
		differences[0]=diffs[0];
	}else{
		*min=new_lim[0];
		*max=new_lim[1];
		differences[0]=diffs[0];
		differences[1]=diffs[1];
	}

	return 0;
}

//double sigma_func(double val,double **par){
//	double x=*(*(par));
//	double Q2=*(*(par)+1);
//	double* sigpar=*(par+1);
//	double* sudpar=*(par+2);
//
//	return(SIGMA(val,x,Q2,sigpar,sudpar));
//}

double sigma_func(double r, double** par){
	return(SIGMA(r,*(*(par)),*(*(par)+1),*(par+1),*(par+2))/(**(par+1))  );
}

//double generate_point(double *par,double Q2,double x,char* filename){
double generate_points(double x, double **param ){
	//double ** param;
	//double *sigpar = *(par+1);
	//double *sudpar = *(par+2);

	double tolerance =1.0e-6;
	double min=1.0e-3;
	double max=9.0;
	double diffs[2];

	*(*(param))=x;

	double point;
	for(int i=0;i<10;i++){
		pick_interval(&sigma_func,param,&min,&max,diffs,0.63212 ); 
		if(fabs(diffs[0])<tolerance|| fabs( diffs[1])<tolerance){
			if(fabs(diffs[0])<fabs(diffs[1]) ){
				point=min;
				continue;
			}else{
				point=max;
				break;
			}
		}
	}
	return 4.0/pow(point,2);
}	


int main(int argc, char ** argv){
	int xlen=100;
	double xarr[xlen+1];
	char file_name[500];
	double param[10];
	double sigpar[10];
	double sudpar[10];
	double x, Q2;
	
	double *par[3];

	read_options(argc,argv,param,&x,&Q2, file_name);
	parameter(param,sigpar,sudpar);
	double var[2];
	var[1]=Q2;
	//*par={x,Q2};
	*(par)=var;
	*(par+1)=sigpar;
	*(par+2)=sudpar;	
	float dum;
	char* end;
	
#if (MODEL==1||MODEL==3)	
	approx_xg(sigpar+1);//generate chebyshev coefficients
#endif
	//printf("gluon ready\n");
	double point=0;
	for (int i=0;i<=xlen;i++){
		*(xarr+i)=pow(10.0,-6+4*(((double)i)/xlen) );
		
	}
	FILE* file=fopen(file_name,"w");
	if(file==NULL){
		printf(" critical line :: file error. %s.\n", file_name);
	}	
	plot(&generate_points,xarr,xlen+1,par,file);
	fclose(file);
	
	return 0;
}



