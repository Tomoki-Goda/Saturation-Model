#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#include"control_tmp.h"
#include"./control-default.h"
#include"./constants.h"

extern double SIGMA(double ,double, double ,double*,double *);
extern void approx_xg(double*);
		
#if MODEL==0
#define SIGMAP sigma_gbw
#elif MODEL==1
#define SIGMAP sigma_bgk
#elif MODEL==2
#define SIGMAP siga_gbs
#else
#define SIGMAP sigma_s
#endif


int pick_interval(double(*func)(double , double **), double ** par , double *min ,double *max,double *differences, double aim){
	if (max<min){
		printf("max %f less than min%f \n", *max,*min);
		return 1;
	}
	int points_n=11;
	double val[points_n];
	double step = (*max-*min)/(points_n-1);
	double diff[points_n];

	double new_lim[2];
	//double diffs[2];
	double *diffs;

	for(int i= 0 ;i<points_n;i++){
		
		val[i]=(*func)(*min+i*step , par);
		//printf("\nsigma = %.3e \t at \t %.3e\n",val[i], *min+i*step) ;		
		diff[i]=(val[i]-aim);
		//printf("diff %d \t %.2e\t",i,diff[i]);
		//printf("prod\t%.3e\n",(diff[i-1]*diff[i]));
		if(i!=0 && (diff[i-1]*diff[i])<0){
			//sign of difference changed -> went pass the point.
			new_lim[0]=*min+(i-1) *step;
			new_lim[1]=*min+i*step;
			//printf("found %.2e %.2e\n",diff[i-1],diff[i]);
			diffs=(diff+i-1);
			break;
		}
		//printf("%d\t%d\n",(i-1),points_n);
		if((i+1)==points_n){
			printf("not found\n");
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

double sigma_func(double val,double **par){
	double x=*(*(par));
	double Q2=*(*(par)+1);
	double* sigpar=*(par+1);
	double* sudpar=*(par+2);

	return(SIGMA(val,x,Q2,sigpar,sudpar));
}


double generate_point(double *par,double Q2,double x,char* filename){
	//FILE* file=fopen(filename, "w");

	double* sigpar= par;
	double *sudpar;
#if MODEL==22||MODEL==2
	sudpar=(par+3);
#elif MODEL==3
	double sudpar[10];
	sudpar[0]=par[3]*par[5] ;//C*C2
	sudpar[1]=par[4]*sqrt(par[5]);//rmax mu02=C/rmax^2
#if SUDAKOV==2
	sudpar[2]=par[6];
	sudpar[3]=par[7];
#endif
#endif
	double tolerance =1.0e-5;

	//double Q2 = 10.0;

	double min=1.0e-1;
	double max=1.0e+1;
	double *param[3];
	double var[2];
	var[0]=x;
	var[1]=Q2;
	*(param)=var;
	*(param+1)=sigpar;
	*(param+2)=sudpar;

	double diffs[2];

	double point;
	for(int i=0;i<100;i++){
		//printf("%.3e \t%.3e\t from %f \n",diffs[0],diffs[1],0.63212*sigpar[0]);
		pick_interval(&sigma_func,param,&min,&max,diffs,0.63212*sigpar[0]); 
		
		if(fabs(diffs[0])<tolerance|| fabs( diffs[1])<tolerance){
			//printf("%.3e\t%.3e\n",diffs[0],diffs[1]);
			if(fabs(diffs[0])<fabs(diffs[1]) ){
				point=min;
				continue;
			}else{
				point=max;
				break;
			}
			
		}
	}
	//printf("point %.5e\n",point);
	return point;
}	


int main(int argc, char ** argv){
	char parfilename[100]="";
	char resultfilename[100]="";
	char name[20]="";
	float par;
	double param[10];//just any number large enough
	float dum;
	//char dumc[100];
	double Q2=100;
	double x=0.001;
	char* end;
	
	for(int i=1;i<=argc/2;i++){
		if( strcmp(argv[2*i-1],"-Q2")==0 ){
			Q2=strtod(argv[2*i],&end);
		}else if( strcmp(argv[2*i-1],"-out")==0 ){
			sprintf(resultfilename,"%s",argv[2*i]);
		}else if( strcmp(argv[2*i-1],"-in" ) ==0){
			sprintf(parfilename,"%s" ,argv[2*i]);
		}else{
			printf("Please Use flag, -Q2, -x, -out , -in\n\n");
		}
			
	}
	printf("Q2= %f \n",Q2);
	
	
	printf("%s\n%s\n",parfilename,resultfilename);
	FILE* parfile=fopen(parfilename,"r");
	
	if(parfile==NULL){
		printf(" file error\n");
		return 1;
	}
	
	fscanf(parfile,"%s\t%f",name,&dum);//line for Qup
	
	for(int i=0;i<N_PAR;i++){
		fscanf(parfile,"%s\t%lf\t%f\n",name,param+i,&dum);
		fprintf(stdout,"%s \t\t %.3e  \n",name ,*(param+i));
	}
	fclose(parfile);
	
	FILE* file=fopen(resultfilename,"w");

#if (MODEL==1||MODEL==3)	
	approx_xg(param+1);//generate chebyshev coefficients
#endif
	//printf("gluon ready\n");
	double point=0;
	for (int i=0;i<100;i++){
		x=pow(10.0,-6+4*(((double)i)/100) );
		//printf("x=%f\n",x); 
		point =generate_point(param, Q2,x ,resultfilename);
		fprintf(file,"%f\t%f\n",x,point);
	}
	fclose(file);
}



