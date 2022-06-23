#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include"./constants.h"
#include"./control-default.h"

//#include"./chebyshev-1.h"
#include"./gluon-chebyshev.h"

#include"./Parameters.h"
#include"./dipole-cross-section.h"


int main(){
	double val, x, q2;
	/*
	for(int i=1;i<=100;i++){
		val=change_var_compactify_log(1,100,i);
		printf("%f\t %f",(double)i, val);
		val=change_var_revert_log(1,100,val);
		printf("\t%f", val);
		printf("\n");
	}
	for(int i=1;i<=100;i++){
		val=change_var_compactify(1,100,i);
		printf("%f\t %f",(double)i, val);
		val=change_var_revert(1,100,val);
		printf("\t%f", val);
		printf("\n");
	}
	
	for(int i=1;i<=100;i++){
		val=change_var_compactify(1.0/100,1, 1.0/i);
		printf("%f\t %f",1.0/(double)i, val);
		val=change_var_revert(1.0/100,1,val);
		printf("\t%f", val);
		printf("\n");
	}*/
	
	double var[2];
	FILE *file=fopen("./xg.txt","w");
	if(file==NULL){
	 printf("error");
	 return 1;
	}
	approx_xg(par_start+1);
	
	
	for(int i=0;i<3;i++){
		//printf( "%f\t%d\n",pow(10,-2) ,i);
		x=pow(10.0,-5+ i);
		//printf( "%f\t%d\n",x,i);
		
		for(int j=0;j<1000;j++){
			q2=j+0.5;
			val=xg_chebyshev(x,q2);
			fprintf(file, "%f\t%f\t%f\t", q2,x,val);
			(*var)=change_var_compactify_log(1.0e-7,1.0,x);
			(*(var+1))=change_var_compactify_log(0.5,1.0e+10,q2);
			
			val=eval_xg(var ,par_start+1);
			fprintf(file, "%f\t%f\t%f\n", val,*var,*(var+1));
			//printf("%f->%f\t ->%f\n",q2,*(var+1),change_var_revert_log(0.05,1.0e+10, *(var+1) ));
		}
	}
	fclose(file);
	return 0;
}





