#include<stdio.h>
#include<stdlib.h>
#include<string.h>

void system_printf(const char* str){
	printf("%s\n",str);
	system(str);
}

int main(int argc, char** argv){
	if(argc<1){
		printf("specify directory");
		return 1;
	}
	char command[1000];
	char dir[50]="";
	char runfile[100]="";
/////////////////   compile   //////////////////////////////////////
	for(unsigned i=1;i<argc;i++){
		strcpy(dir,argv[i]);
		sprintf(command,"rm %s/main",dir);
		system(command);
		
		printf("going to %s/\n",dir);
		sprintf(command, "cp \"%s/control.h\" ./control_tmp.h",dir);
		system(command);
		
		sprintf(command, "tar -czvf %s/control-default.h.tar.gz control-default.h ",dir);
		system(command);
		
		//system("make");
		//sprintf(command,"./main \"%s/result.txt\"",dir);
		
		sprintf(command,"gcc main.c -o %s/main -lm -lmathlib -lkernlib -lpacklib",dir);
		system(command);
		
		
		sprintf(command,"rm \"./control_tmp.h\"");
		system(command);
		
		
	}
	
///////////////////////   prepare command   /////////////////////
	if(argc>2){
		//strcpy(command,"parallel -keep-order --link -j 2 --line-buffer :::");
		strcpy(command,"parallel --link -j 2 --line-buffer :::");
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
	if(argc>2){
		strcat(command, " :::");
	}
	
	for(unsigned i=1;i<argc;i++){
		strcpy(dir,argv[i]);
		sprintf(runfile," %s",dir);
		strcat(command, runfile);
	}
	
	/*
	if(argc>2){
		strcat(command, " :::");
	}
	
	for(unsigned i=1;i<argc;i++){
		strcpy(dir,argv[i]);
		sprintf(runfile," %s/log.txt",dir);
		strcat(command, runfile);
	}
	*/
////////////////   RUN!   //////////////////////////
	system_printf(command);
	
	return 0;

}
