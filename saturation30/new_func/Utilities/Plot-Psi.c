#include<math.h>
#include<stdio.h>
#include"control.h"
#include"../control-default.h"
#include"../constants.h"
//#include"./photon-wave-function-2.h"
#include<stdlib.h>
#include<unistd.h>

extern double psisq_z_int(double ,double ,int);

int main(int argc, char* argv[]){
	FILE *file;
	
	unsigned n=100;
	double val;
	double r;
	double Q2=0.0;
	//char name[100];
	
	int opt=0;
	for(int i =1;i<argc;i++){
		printf("%s ", argv[i]);
	}
	printf("\n");

	while((opt=getopt(argc,argv,"Q:o:") )!=-1 ){
		//opt=getopt(argc,argv,":Q2:out:");
		printf("%d", opt);
		switch(opt ){
			case 'o':
			//	break;
			//case 'u':
			//	break;
			//case 't':
				printf("%s\n", optarg);
				file=fopen(optarg,"w");
				break;
			case 'Q':
			//	break;
			//case '2':
				Q2=atoi(optarg);
				printf("%.2e\n",Q2);	
				break;
			//case 'x':
			//	x=atod(optarg);
			//	break;
			case '?':
				printf("options should be -Q2 , -x , -h\n");
				break;

		}
	}
	//sprintf(name,"./photon-psi%s.txt",argv[1]);
	
	//file=fopen(name,"w");
	for (unsigned i=0 ; i<=n;i++){
		r=5*pow(10,-2+3*((double)i)/n);
		val=0;
		for(int j =0 ;j<(NF-1);j++){
			val+=pow(r,-2)*psisq_z_int(r,Q2,j);
		//	val+=psisq_z_int(r,Q2,j);
		}
		//fprintf(file,"%f\t%f\n",pow(r,-2),val);
		fprintf(file,"%f\t%f\n",r,val);
	}
	fclose(file);
}
