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
		sprintf(command,"rm %s/Plot",dir);
		system(command);
		
		printf("going to %s/\n",dir);
		sprintf(command, "cp \"%s/control.h\" ./control_tmp.h",dir);
		system(command);
		
		sprintf(command,"gcc Plot-DP.c -o %s/Plot -lm -lmathlib -lkernlib -lpacklib",dir);
		system(command);
		
		sprintf(command,"rm \"./control_tmp.h\"");
		system(command);
		
		sprintf(command,"%s/Plot %s",dir,dir);
		system(command);
	}
	
	
	return 0;

}
