#include<stdio.h>
#include<stdlib.h>
///syntax///
// ./main <no. of parameters > <input files...> <output file>
//#define line_n 3

int main(int argc, char** argv){
	//printf("%s",argv[argc-1]);
	FILE* resfile;
	FILE* outfile=fopen(argv[argc-1],"w");
	
	
	char name[20];
	float value;
	float error;
	float chisq;
	int ndata;
	
	//int line_n= *(argv[1])-'0';//char to int
	int line_n= atoi(argv[1]);
	printf("%d parameters \n",line_n);
		
	for(unsigned i =2 ;i<argc-1;i++){	
		resfile=fopen(argv[i],"r");
		//fprintf(outfile,"%s",argv[i]);
		
		fscanf(resfile,"%s\t%f",name, &value);
		fprintf(outfile,"%.0f",value);
		for(unsigned line =0; (line<(line_n)) ; line++){
			fscanf(resfile,"%s\t%f\t%f\n",name,&value, &error );
			fprintf(outfile,"& %.2e {\\tiny $\\pm$ %.2e }",value, error);
			//fprintf(stdout,"& %f $\\pm$ %f \n",value, error);
		}
		fscanf(resfile,"%s\t%f\t%f\n",name,&chisq, &error );
		fprintf(outfile,"& %.2e  ",chisq);
		fscanf(resfile,"%s\t%d\n",name,&ndata);
		fprintf(outfile,"/ %d  ",ndata);
		
		fprintf(outfile,"= %.2f  ",chisq/ndata);
		fprintf(outfile,"\\\\ \\hline \n");
			
		fclose(resfile);	
	}
	fclose(outfile);
}
