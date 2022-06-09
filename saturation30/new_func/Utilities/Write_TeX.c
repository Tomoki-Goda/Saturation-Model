#include<stdio.h>
#include<stdlib.h>
#include<string.h>

///syntax///
// ./Write_TeX <input files...> <output file>
//#define line_n 3

int main(int argc, char** argv){
	//printf("%s",argv[argc-1]);
	FILE* resfile;
	FILE* outfile=fopen(argv[argc-1],"w");
	
	
	char name[20];
	float value;
	float error;
	//float chisq;
	int ndata;
	unsigned count=0;
	
	//int line_n= *(argv[1])-'0';//char to int
	//int line_n= atoi(argv[1]);
	
	//printf("%d parameters \n",line_n);
	
	
	for(unsigned i =1 ;i<argc-1;i++){	
		count=0;
		
		resfile=fopen(argv[i],"r");
		//fprintf(outfile,"%s",argv[i]);
		if(resfile==NULL){ 
			printf("makelist error no file");
			return 1;
		}
		fscanf(resfile,"%s\t%f",name, &value);//Qup
		fprintf(outfile,"%.0f",value);
		
		
		//for(unsigned line =0; (line<(line_n)) ; line++){
		for(unsigned line =0; line<10 /* > max possible no. of parameter*/; line++){
			fscanf(resfile,"%s\t%f\t%f\n",name,&value, &error );
			if(strcmp(name,"chisq")==0){
				printf("%d parameters \n",count);
				break;	
			}
			fprintf(outfile,"& %.2e {\\tiny $\\pm$ %.2e }",value, error);
			count++;
			//fprintf(stdout,"& %f $\\pm$ %f \n",value, error);
		}
		//fscanf(resfile,"%s\t%f\t%f\n",name,&chisq, &error );
		//fprintf(outfile,"& %.2e  ",chisq);
		fprintf(outfile,"& %.2e  ",value);
		fscanf(resfile,"%s\t%d\n",name,&ndata);
		fprintf(outfile,"/ %d  ",ndata);
		
		fprintf(outfile,"= %.2f  ",value/(ndata-count));
		fprintf(outfile,"\\\\ \\hline \n");
			
		fclose(resfile);	
	}
	fclose(outfile);
	return 0;
}
