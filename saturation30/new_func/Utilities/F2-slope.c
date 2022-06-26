#include"./f2.h"


double f2_slope(double Q2, double**pars){
	double  err;
	*(*(pars)+1)=Q2;
	//printf("x== %.2e\t Q2== %.2e\n", pars[0][0],pars[0][1]);
	//getchar();
	double xdiff=0.01;
	double logx=log10(**pars);
	double f2[2];
	simpson1dA(&f2_integrand,pars,1.0e-5,30,100,f2,&err);
	(*(*pars)) = pow(10.0, logx+xdiff);
	simpson1dA(&f2_integrand,pars,1.0e-5,30,100,f2+1,&err);
	double slope = (log10(f2[0]) - log10( f2[1]) )/xdiff;
	//printf("%f\n",res);
	(*(*pars)) = pow(10.0, logx);
	return slope;
}



int main (int argc, char** argv){
	char file_name[500];
	int Q2len=1000;
	double Q2arr[Q2len+1];
	double param[10];
	double x=0;
	double Q2=0;
	
	read_options(argc,argv,param,&x,&Q2,file_name );
	//printf("x= %.2e\tQ2= %.2e\n",x,Q2);
	
	for(int i =0;i<=Q2len;i++){
		*(Q2arr+i)=pow(10.0, -1+3*((double)(i))/Q2len);
	}

	double *par[3];
	double var[2];
	var[0]=x;
	var[1]=0;
	
	//printf("x= %.2e\tQ2= %.2e\n",var[0], var[1]);
	//getchar();
	par[0]=var;
	double sigpar[10];
	double sudpar[10];
	parameter(param,sigpar,sudpar);
	par[1]=sigpar;
	par[2]=sudpar;

#if (MODEL==1||MODEL==3)	
	approx_xg(sigpar+1);//generate chebyshev coefficients
#endif

	FILE *file=fopen(file_name,"w");
	if(file==NULL){
		printf("F2-slope:: file error. %s.\n", file_name);
	}	
	plot(&f2_slope,Q2arr,Q2len+1,  par,  file);
	fclose(file);

}



