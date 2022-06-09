#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

//Syntax is ./Auto-Plot-F2 <log10 x> < directories of result.txt & control.h>

void system_printf(const char* str){
	printf("%s\n",str);
	system(str);
}

int main(int argc, char** argv){
	if(strcmp(argv[1],"-h")==0||strcmp(argv[1],"--help")==0){
		printf( "Syntax is ./Auto-Plot-F2 <-log10 x> < directories of result.txt & control.h>\n");
		printf( "e.g. \n\t ./Auto-Plot-F2 3 ./Directory/* \nFor multiple folders in Directory\n");
		return 0;
	}
	
	if(argc<3){
		printf("specify directory");
		return 1;
	}
	char command[1000];
	char dir[50]="";
	char runfile[100]="";
	
	//int Q2=100;
	int xpow=3;
	char *end;
	//Q2=strtod(argv[1],&end);
	xpow=strtod(argv[1],&end);
/////////////////   compile   //////////////////////////////////////
	for(unsigned i=2;i<argc;i++){
		strcpy(dir,argv[i]);
		sprintf(command,"rm %s/Plot",dir);
		system(command);
		
		printf("going to %s/\n",dir);
		sprintf(command, "cp \"%s/control.h\" ./control_tmp.h",dir);
		system(command);
		
		sprintf(command,"gcc ./Utilities/Plot-F2.c -o %s/Plot-F2 -lm -lmathlib -lkernlib -lpacklib",dir);
		system(command);
		
		sprintf(command,"rm \"./control_tmp.h\"");
		system(command);
		
		printf("%s/Plot-F2 -in \"%s/result.txt\"  -out \"%s/pointsF2-x%d.txt\"   -x %f\n",dir,dir,dir,xpow, pow(10.0,-xpow) );
		
		sprintf(command,"%s/Plot-F2 -in \"%s/result.txt\"  -out \"%s/pointsF2-x%d.txt\"  -x %f\n",dir,dir,dir,xpow,pow(10.0,-xpow) );
		
		system(command);
	}
	
	
	return 0;

}