#include"./f2.h"


double f2(double Q2, double**pars){
	double res=0, err=0,val=0;
	*(*(pars)+1)=Q2;
	//printf("x== %.2e\t Q2== %.2e\n", pars[0][0],pars[0][1]);
	for(int i =0;i<(NF-1);i++){
		pars[3][0]=i;
		res=0;
		simpson1dA(&f2_integrand,pars,1.0e-5,30,250,&res,&err);
		val+=res;
	}
	//printf("%f\n",res);
	return val;
}



int main (int argc, char** argv){
	char file_name[500];
	int Q2len=30;
	double Q2arr[Q2len+1];
	double param[10];
	double x=0;
	double Q2=0;
	
	read_options(argc,argv,param,&x,&Q2,file_name );
	printf("x= %.2e\tQ2= %.2e\n",x,Q2);

	for(int i =0;i<=Q2len;i++){
		*(Q2arr+i)=pow(10.0, -1+3*((double)(i))/Q2len);
	}

	double *par[4];
	double var[2];
	var[0]=x;
	var[1]=0;
	
	//printf("x= %.2e\tQ2= %.2e\n",var[0], var[1]);
	par[0]=var;
	double sigpar[10];
	double sudpar[10];
	parameter(param,sigpar,sudpar);
	par[1]=sigpar;
	par[2]=sudpar;
	double flavour=0;
	par[3]=&flavour;
#if (MODEL==1||MODEL==3)	
	approx_xg(sigpar+1);//generate chebyshev coefficients
#endif

	FILE *file=fopen(file_name,"w");
	if(file==NULL){
		printf("F2:: file error. %s.\n", file_name);
	}	
	plot(&f2,Q2arr,Q2len+1,  par,  file);
	fclose(file);

}



