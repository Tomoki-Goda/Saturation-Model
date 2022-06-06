#include<stdio.h>


int main(int argc, char** argv){
	char addline[100];
	sprintf(addline ,"%s" ,argv[1]);
	printf("%s\n" ,addline);
	
	FILE* file;
	char filename[100];
	for(unsigned i=2;i<argc;i++){
		sprintf(filename, "%s/control.h",argv[i]);
		file=fopen(filename,"a");
		if(file==NULL){
		 	printf("file error");
		 	return 1;
		 }
		 fprintf(file,"%s\n",addline);
		 fclose(file);
	}

}
