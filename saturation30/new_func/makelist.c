#include<stdio.h>


int main(int argc, char** argc){
	FILE* readfile;
	char name[argc][20];
	double value[argc];
	double error[argc]
	for(unsigned i=1;i<argc;i++){
		readfile=fopen(argc[i],"r");
		while(true){
		
			fscanf(readfile,"%s\t%f\t%f\n",(name+i) ,vale,error)
	}

}
