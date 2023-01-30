#include<stdio.h>
#include<unistd.h>
#include<getopt.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>

typedef struct{
	char input_file_name[100];
	char output_file_name[100];
	double Q2;
	double x;
	double k;
	int fl;
	double W;
	char data_file_name[100];
	double xmax,xmin,beta;
} options;

options read_options(int argc, char** argv  ){
	int c, indexptr;
	options data;
	struct option long_options[]={
		{"Q2",required_argument,0,'q'},
		{"k",required_argument,0,'k'},
		{"x",required_argument,0,'x'},
		{"flavour",required_argument,0,'f'},
		{"in",required_argument,0,'i'},
		{"out",required_argument,0,'o'},
		{"W",required_argument,0,'w'},
		{"data",required_argument,0,'d'},
		{"xmax",required_argument,0,'M'},
		{"xmin",required_argument,0,'m'},
		{"beta",required_argument,0,'b'}

	};
	char shortopts[]="q:k:x:f:i:o:w:d:M:m:b:";
	
	while(1){
		c=getopt_long_only(argc, argv,shortopts,long_options,&indexptr );
		if(c==-1){break;}

		switch(c){
			case 'q':
				data.Q2=atof(optarg);
				break;
			case 'k':
				data.k=atof(optarg);
				break;
			case 'x':
				data.x=pow(10,-atoi(optarg));
				break;
			case 'f':
				data.fl=atoi(optarg);
				break;
			case 'i':
				strcpy(data.input_file_name,optarg);
				break;
			case 'o':
				strcpy(data.output_file_name,optarg);
				break;
			case 'w':
				data.W=atof(optarg);
				break;
			case 'd':
				strcpy(data.data_file_name,optarg);
				break;
			case 'M':
				data.xmax=atof(optarg);
				break;
			case 'm':
				data.xmin=atof(optarg);
				break;
			case 'b':
				data.beta=atof(optarg);
				break;
			default:
				printf("Unknown option\n");
		};

	}
	return(data);
}

//int read_parameters(FILE * parfile, double *param){
int read_parameters(FILE * parfile, std::vector< double>& param){
	double dum;
	char name[20];
/////////////////////param result data //////////////////////////////
	fscanf(parfile,"%s\t%lf",name,&dum);//line for Qup
	
	//for(unsigned i=0;i<N_PAR;i++){
	for(int i=0;i<10/*any large enough*/;i++){
		fscanf(parfile,"%s\t%lf\t%lf\n",name,&(param[i]),&dum);
		if(strcmp(name,"chisq")==0){
			break;
		}
		fprintf(stdout,"%s \t\t %.3e  \n",name ,(param[i]));
	}
	
	return 0;

}

//int main(int argc , char**argv){
//	struct info data;//structure to store data;
//	read_options(argc,argv,&data);
//	printf("%f %f\n",data.Q2, data.x);
//
//	return 0;
//}

