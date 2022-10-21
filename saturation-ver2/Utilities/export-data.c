
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<cfortran.h>
#include<minuit.h>

#include"control.h"
#include"control-default.h"
#include"constants.h"



//#include"./Parameters.h"
//#include"./f2.h"
//#define MAXN 600

extern double SIGMA_PREC;
static unsigned N_DATA;
extern int N_SIMPS;
extern int N_CHEB;


extern int load_data(void);
extern void generate_psi_set(void);
FILE* log_file;

extern void approx_xg(const double *);
extern int parameter(const double*,double *, double*);

#include"../Utilities/plot.c"

//extern void fcn(const int *npar, const double grad[], double*fcnval, const double *par,const unsigned *iflag,void (*dum)(void) );
//extern void dum_func(void);
//extern void save_f2(char* file_name);
extern void  export_data(FILE * file,double *sigpar,double *sudpar);

void log_printf(FILE* file,char* line){
	if(file!=stdout){
		fprintf(file,"%s",line);
	}
#if (PRINT_PROGRESS==1)
	fprintf(stdout,"%s",line);
#endif
}

int main(int argc, char** argv){
	char file_name[500];
	double param[10]={0};
	double sigpar[10];
	double sudpar[10];
	double x, Q2;
	log_file=stdout;
	
	
	read_options(argc,argv,param,&x,&Q2, file_name);
	

	//static double sigpar[10],sudpar[10];	
	parameter(param,sigpar,sudpar);
#if (MODEL==1||MODEL==3)
	approx_xg(sigpar+1);//generate chebyshev coefficients
#endif
	
	N_DATA=load_data(); 
	
	FILE* file=fopen(file_name,"w");
	export_data(file, sigpar, sudpar);
	fclose(file);
	
	
	return 0;
	

}
