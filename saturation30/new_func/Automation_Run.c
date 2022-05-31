#include<stdio.h>
#include<stdlib.h>
#include<string.h>


int main(int argc, char** argv){
	if(argc<1){
		printf("specify directory");
		return 1;
	}
	char command[1000];
	char dir[50]="";
	char runfile[100]="";

	for(unsigned i=1;i<argc;i++){
		strcpy(dir,argv[i]);
		
		printf("going to %s/\n",dir);
		sprintf(command, "cp \"%s/control.h\" ./control_tmp.h",dir);
		system(command);
		//system("make");
		//sprintf(command,"./main \"%s/result.txt\"",dir);
		
		sprintf(command,"gcc main.c -o %s/main -lm -lmathlib -lkernlib -lpacklib",dir);
		system(command);
		
		
		sprintf(command,"rm \"./control_tmp.h\"");
		system(command);
	
		
	}
	
	
	if(argc>2){
		strcpy(command,"parallel --link -j 4 :::");
	}else{
		strcpy(command,"");
	}
	for(unsigned i=1;i<argc;i++){
		strcpy(dir,argv[i]);
		sprintf(runfile," %s/main",dir);
		strcat(command, runfile);
		//sprintf(command,"./%s/main \"%s/result.txt\"&",dir,dir);
		//system(command);
	}
	strcat(command, " :::");
	for(unsigned i=1;i<argc;i++){
		strcpy(dir,argv[i]);
		sprintf(runfile," %s/result.txt",dir);
		strcat(command, runfile);
		//sprintf(command,"./%s/main \"%s/result.txt\"&",dir,dir);
		//system(command);
	}
	printf("\n\n %s \n\n",command);
	system(command);
	
	return 0;

}
