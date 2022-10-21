#include"./ww-gluon.h"

int main (int argc, char** argv){

	double val;
	FILE *file;
	char file_name[500];
	//double k2, x , q2;
	double param[10];
	double sudpar[10];
	double sigpar[10];
	//double step=((double)R_MAX)/(2*n);
	double K;
	
	
	read_options(argc,argv,param,&(ww_par.X),&(ww_par.Q2), file_name);
	parameter(param,sigpar,sudpar);
	ww_par.SIGPAR=sigpar;
	ww_par.SUDPAR=sudpar;

	printf("%.3e %.3e\n",ww_par.X, ww_par.Q2);
		
#if (MODEL==1||MODEL==3)	
	approx_xg(sigpar+1);//generate chebyshev coefficients
#endif
	file=fopen(file_name,"w");
	if(file==NULL){
			printf("tmd-gluon:: file can't be opened. %s\n",file_name);
			return 1;
	}
	
	for(int i=0;i<25; i++){
		K=pow(10,-1+2*((double)i)/25);
		ww_par.K=K;
		val=ww_integral();
		//val=ww_grad();
		//val*=3.0/(4*PI);
		val*=4.0/(6.0* PI*PI*PI*PI);		
		fprintf(file,"%.5e\t%.5e\n",K*K,val);
	}
	fclose(file);
	
	return 0;
	
}
