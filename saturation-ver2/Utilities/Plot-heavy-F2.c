#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"control.h"
#include"control-default.h"
#include"constants.h"

#include"./options.h"


//#include"./kahnsum.h"
/////////////////////kahn.h/////////////////////
extern double KBN_sum(const double *arr,int len);
extern double kahn_sum(const double *arr,int len);
////////////////chebyshev-1.h/////////////////////////////
double cheb_c(const double * sample_arr, const unsigned* ind1,const unsigned *degree, unsigned dim );
double change_var_revert(double min,double max, double val);
double change_var_revert_log(double min,double max, double val);
double change_var_revert_frac(double , double, double val);
void sample(double func(const double *,const double* ), const double * par,  const unsigned * degree, unsigned dim,  double* sample_arr);


///////////////photon-wave-function.c/////////////////
extern double psisq_z_int(double, double ,int);
extern double mod_x(double,double, int);
//////////////gluon-chebyshev.h///////////////////////
extern void approx_xg(const double *);
/////////////dipole-cross-section.c/////////////////
extern int parameter(const double*,double*,double*);
extern double SIGMA(double , double ,double , const double *, const double*);
///////////Set in main.c//////////////////
extern void log_printf(FILE*,char*);
extern FILE* log_file;

#include"./f2.h"

double DATA[100][4];

struct info data;

int main(int argc , char** argv){
	read_options(argc,argv,&data);
	char input_file_name[100];
	strcpy(input_file_name,data.input_file_name);
	char output_file_name[100];
	strcpy(output_file_name,data.output_file_name);
	int fl=data.fl;
	PLOT_FLAVOUR=fl;
	double pars[10];
	char DATA_DIR[100];
	if(fl==2){
		strcpy(DATA_DIR,"./data/H1F2c.txt");
	}else if(fl==3){
		strcpy(DATA_DIR,"./data/H1F2b.txt");
	}else{
		printf("choose flavour 2 or 3\n");
	}
	
	FILE* input_file=fopen(input_file_name,"r");
	read_parameters(input_file,pars);
	static double sigpar[10],sudpar[10];	
	parameter(pars,sigpar,sudpar);
#if (MODEL==1||MODEL==3)
	approx_xg(sigpar+1);//generate chebyshev coefficients
#endif
	
	fclose(input_file);	
//////////////////////////////////////////////
/// A bit problematic part, use the name for output.
/// insert "-data" before ".txt" for the datafile.
////////////////////////////////////////////////	
	char data_file_name[100];
	strcpy(data_file_name,output_file_name);
	printf("out= %s\n",output_file_name);
	char newend[]="-data.txt";
	int pos, endlen=strlen(newend);
	printf("%d\n", endlen);
	pos=strlen(data_file_name);
	for(int i=0;i<endlen+1;i++){
		data_file_name[pos-4+i]=newend[i];		
	}
	//strcat(data_file_name+((int)strlen(output_file_name)-4), "-data.txt" );
	printf("%s\n",data_file_name);
	
////////////////////////////////////////////////////	
	
	FILE* data_file = fopen(data_file_name, "w");
	
	FILE* output_file = fopen(output_file_name, "w");
	
	
	double dummy,x,q2,val;
	char name[3];
	FILE* file=fopen(DATA_DIR,"r");
	fscanf(file, "%*[^\n]");
	int j=-1;
	while(!feof(file)&& ++j<100 ){
		fscanf(file ,"%s\t%lf\t%lf\t%lf\t%lf\t%lf\t",name,&(DATA[j][0]),&(DATA[j][1]),&dummy,&dummy,&(DATA[j][2]) );
		fscanf(file ,"%lf\t%lf\t%lf\t",&dummy,&dummy,&(DATA[j][3]) );
		for(int i=0; i<11;i++){
			fscanf(file,"%lf",&dummy);
		}
		DATA[j][3]*=(DATA[j][2]/100);
		//printf("\t data %d %s Q2=%.5e,x=%.5e, F2=%.5e +- %.5e \n" , j, name,DATA[j][0],DATA[j][1],DATA[j][2],DATA[j][3]);
		if(feof(file)){
			break;
		}
		if(DATA[j][1]>1.0e-2){
			continue;
		}
		fprintf(data_file,"%.5e\t%.5e\t%.5e\t%.5e\n",DATA[j][0],DATA[j][1],DATA[j][2],DATA[j][3]);
		printf("data %d %s Q2=%.5e,x=%.5e, F2=%.5e +- %.5e \n" , j, name,DATA[j][0],DATA[j][1],DATA[j][2],DATA[j][3]);
		for(int i=0;i<5;i++){
			x=(0.5+(3.5)*((double)i)/5)*DATA[j][1];
			q2=DATA[j][0];
			//printf("%e %e %e %e\n",x, q2,sigpar[0],sudpar[0]);
			val=f2_2(x,q2, sigpar , sudpar);
			fprintf(output_file,"%.5e\t%.5e\t%.5e\n",DATA[j][0],x,val);
		}
	}
	fclose(output_file);
	fclose(data_file);		
	return 0;
}
