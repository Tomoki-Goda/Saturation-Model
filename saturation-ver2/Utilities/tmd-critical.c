#include"./tmd-gluon.h"

double saturation(double step){
	double k_step=0.1;
	double k_min=0.5;
	double k=k_min;
	double prev=1;
	double val;
	
	for(int i=0;i<6;i++){
		//printf("%.3e\t%.3e\n",k_min,k_step );
		k=k_min;
		
		for(int j=0;j<100;j++){
			
			val=grad_k(k,step);
			//printf("%d : %.3e, %.3e, %.3e, %.3e\n",j ,k,k_min, val, prev );
			
			if(j!=0 ){
				if(prev*val<0){
					k_min=k-k_step;
					//printf("\n%d : %.3e, %.3e, %.3e, %.3e\n",i ,k,k_min, val, prev );
					break;
				}
			}
			//printf("%f, %f, %f\n",k, val, prev );
			prev=val;
			k+=k_step;
		}
		if(i!=5){
			k_step*=0.05;
		}
		
	}
	//printf("\nk=%.3e, val= %.3e, prev= %.3e\n\n",k, val, prev );
	return(k - k_step/2 );
}

int main (int argc, char** argv){

	double val;
	FILE *file;
	char file_name[500];
	double k, x , Q2;
	double param[10];
	double sudpar[10];
	double sigpar[10];
	double step=(60.0)/(2*n);
	
	read_options(argc,argv,param,&x,&Q2, file_name);
	parameter(param,sigpar,sudpar);
	printf("%.3e %.3e\n",x, Q2);
	
#if (MODEL==1||MODEL==3)	
	approx_xg(sigpar+1);//generate chebyshev coefficients
#endif

	
	file=fopen(file_name,"w");

	if(file==NULL){
		printf("tmd-gluon:: file can't be opened. %s\n",file_name);
		return 1;
	}
	for (int i=0; i<50; i++){
		x= pow(10,-5+((double)3*i)/50);
		sample_sigma( sample ,  step,  x, Q2, sigpar,  sudpar);
		
		val= saturation(step);
		//val*=k*k;
		printf("%.5e\t%.5e\t%.5e\n",x, val, grad_k(val,step));
		
		fprintf(file,"%.5e\t%.5e\n",x, val*val);
	}
	fclose(file);
	
	return 0;
}

