#include"./ww-gluon.h"

int main (int argc, char** argv){

	double val;
	FILE *file;
	char file_name[500];
	//double k2, x , q2;
	//double x;
	double param[10];
	//double sudpar[10];
	//double sigpar[10];
	//double step=((double)R_MAX)/(2*n);
	
	
	
	read_options(argc,argv,param,&X,&Q2, file_name);
	parameter(param,SIGPAR,SUDPAR);
	printf("%.3e %.3e\n",X, Q2);
		
#if (MODEL==1||MODEL==3)	
	approx_xg(SIGPAR+1);//generate chebyshev coefficients
#endif
	file=fopen(file_name,"w");
	if(file==NULL){
			printf("tmd-gluon:: file can't be opened. %s\n",file_name);
			return 1;
	}
	
	for (int i=0; i<=5; i++){
		X= pow(10,-6+((double)4*i)/5);
		
		val= ww_saturation();
		printf("%.5e\t%.5e\t%.5e\n",X,K, 4*val*val);
		fprintf(file,"%.5e\t%.5e\n",X, 4*val*val);
	}
	fclose(file);
	
	return 0;
	
}