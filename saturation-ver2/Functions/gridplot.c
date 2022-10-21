#include<stdio.h>
#include<math.h>
#include<time.h>
#include<cfortran.h>
#include<minuit.h>

#include"control.h"
#include"control-default.h"
#include"constants.h"

//#include"./Parameters.h"
//#include"./main.h"
//#define MAXN 600

extern double SIGMA_PREC;
static unsigned N_DATA;
extern int N_SIMPS;
extern int N_CHEB;


extern int load_data(void);
extern void generate_psi_set(void);
FILE* log_file;


#include"../Utilities/plot.c"

extern void fcn(const int *npar, const double grad[], double*fcnval, const double *par,const unsigned *iflag,void (*dum)(void) );
extern void dum_func(void);
extern void save_f2(char* file_name);

void log_printf(FILE* file,char* line){
	if(file!=stdout){
		fprintf(file,"%s",line);
	}
#if (PRINT_PROGRESS==1)
	fprintf(stdout,"%s",line);
#endif
}
int main(int argc, char** argv){

	int ndata= load_data();
	
	
		
	
}
