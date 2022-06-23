#include<stdio.h>
#include<string.h>


int main(int argc, char** argv){
	for(unsigned i=1;i<argc;i++){
		if(strcmp(argv[i],"-h")==0|| strcmp(argv[i],"--help")==0){
			printf("syntax is :\n\t <line to add> <directories...>\n");
			printf("e.g. \n\t \"#define PRINT_PROGRESS 1\" ./BGK/M* \n");
			return 0;
		}
		
	}
	
	//if(strcmp(argv[1],"--remove")==0){	
	//}
	
	char addline[1000];
	sprintf(addline ,"%s" ,argv[1]);
	printf("%s\n" ,addline);
	
	FILE* file;
	char filename[1000];
	for(unsigned i=2;i<argc;i++){
		sprintf(filename, "%s/control.h",argv[i]);
		file=fopen(filename,"a");
		if(file==NULL){
		 	printf("file error\n");
		 	return 1;
		 }
		 fprintf(file,"%s\n",addline);
		 fclose(file);
	}
	return 0;
}
