#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"./chebyshev-1.h"

double func(double *var,double *par){
	double val=exp(-3*pow((var[2]-var[0])/var[1],2) );
	//printf("return %f %f %f \n", val,var[0],var[1]);
	
	return val;
}

int main(){
	/////////////////test index////////////
	int indices[5]	={0,0,0,0,0};
	int max[5]	={3,3,3,3,3};
	for(unsigned i=0;i<32;i++){
		for(unsigned j=0;j<5;j++){
			printf("%d",*(indices+j));	
		}
		printf("\t");
		printf("%d ?= %d\n",i, convert_index(indices , max, 5));
		ind_vec_increment(indices , max, 5);
	}
	/////////////////////////////////////
	/////////////Test chebyshev//////////
	
	double Ti[8];
	FILE* file;
	file= fopen("./testplotcheb.txt","w");
	if(file==NULL){
		printf("error\n");
		return(1);
	}
	double x;
	for(unsigned i=0;i<=50;i++){
		x=-1+2*((double)i)/50;
		chebyshevT(x,8,Ti);
		fprintf(file,"%f\t",x);
		for(unsigned j=0;j<8;j++){
			fprintf(file,"%f\t",Ti[j]);			
		}
		fprintf(file,"\n");
	}	
	fclose(file);
	///////////////////////////////////
	
	
	double * dum;
	int dim = 3;
	int degree[3]={25,25,25};
	
	//double coeff[10000];
	int len=1;
	for(unsigned i=0;i<dim;i++){
		len*=degree[i];
	}
	double coeff[len];
	
	cheb_coeff(&func , dum ,degree, dim, coeff  );
	
	file;
	file=fopen("./testchebT.txt","w");
	fclose(file);
	
//	FILE* file;
	file= fopen("./testplot.txt","w");
	if(file==NULL){
		printf("error\n");
		//fclose(file);
		return(1);
	}
	//fprintf(file,"hello\n");
//	printf("open\n");
	
	double var[3]={0,1,-0.5};
	double val=0;
	double val2=0;
	
	for(unsigned k=0;k<2;k++){
		(*(var+2))+=0.25;
		for(unsigned j=0;j<2;j++){
			(*(var+1))-=0.2;
			
			for(unsigned i=0;i<500;i++){
	 			(*var)=-1.0+(2.0* ((double)i)/500);
	 			
				val=chebyshev(degree , dim,  coeff , var);
				
				val2=func(var,dum);
	 			
	 			fprintf(file,"%f\t%f\t%f\n",*var,val,val2);
	 			//fprintf(stdout,"%f\t%f\t%f\n",*var,val,val2);
	 		}
	 		
		}
	}
	fclose(file);
	
	return 0;
}





