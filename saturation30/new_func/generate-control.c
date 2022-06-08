#include<stdio.h>
#include<string.h>
#include<stdlib.h>




int main(int argc, char** argv){
	//unsigned i;
	char filenames[2*5*3*3*2][200];
	
	char def[10]="Default";
	
	char* lmass[2];
	unsigned lmasslen=1;
	*lmass=def;
	char* qup[5];
	unsigned quplen=1;
	*qup=def;
	char* sudakov[3];
	unsigned sudakovlen=1;
	*sudakov=def;
	char* model[3];
	unsigned modellen=1;
	*model=def;
	char* rfix[2];
	unsigned rfixlen=1;
	*rfix=def;
	char *dir;
	
	//this function is used to generate control.h for different settings.
	for(unsigned i=1;i<argc-1; i++){
		if(strcmp(argv[i],"-h" )==0||strcmp(argv[i],"--help" )==0){
			printf("********************************************************************************\n");
			printf("Hello there! Asked for help? Here are some tips.\n");
			printf("Possible options are:\n -lmass, -qup, -model, -sudakov\nList the values or whatever they have to be.\n");
			printf("Control file of every combination will be generated\n" ) ;
			printf("Good luck! and cheerio!\n");	
			printf("********************************************************************************\n");
			return 0;		
		}
		if(strcmp(argv[i],"-lmass" )==0){
			for (lmasslen=0;(i< (argc-1))&&(argv[i+1][0]!='-') ;lmasslen++ ){
				*(lmass+lmasslen)= argv[++i];
			}		
		}else if(strcmp(argv[i],"-qup" )==0){
			for (quplen=0;(i< (argc-1))&&(argv[i+1][0]!='-');quplen++ ){
				*(qup+quplen)=argv[++i];						
			}					
		}else if(strcmp(argv[i],"-model" )==0){
			for (modellen=0;(i< (argc-1))&&(argv[i+1][0]!='-');modellen++ ){
				*(model+modellen)=argv[++i];								
			}				
		}else	if(strcmp(argv[i],"-sudakov" )==0){
			for (sudakovlen=0;(i< (argc-1))&&(argv[i+1][0]!='-');sudakovlen++ ){
				*(sudakov+sudakovlen)=argv[++i];									
			}				
		}else if(strcmp(argv[i],"-rfix" )==0){
			for (rfixlen=0;(i< (argc-1))&&(argv[i+1][0]!='-');rfixlen++ ){
				*(rfix+rfixlen)=argv[++i];								
			}
		}else if(strcmp(argv[i],"-dir" )==0){
			dir=argv[++i];
		}else{
			printf("Unknown option?? %s\n",argv[i] );
			return 1;
		}				
	}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//unsigned totallen=(++lmasslen)*(++quplen)*(++sudakovlen)*(++modellen)*(++rfixlen);
	char command[500];
	unsigned totallen=(lmasslen)*(quplen)*(sudakovlen)*(modellen)*(rfixlen);
	unsigned lm, q,s,mod,rf;
	
	FILE* file;
	
	for(unsigned i=0;i<totallen;i++){
		//Now write files...
		lm=i/(quplen*sudakovlen*modellen*rfixlen);
		q= (i-(lm*quplen*sudakovlen*modellen*rfixlen))/( sudakovlen*modellen*rfixlen);
		s= (i-(lm*quplen*sudakovlen*modellen*rfixlen)-(q*sudakovlen*modellen*rfixlen) )/(modellen*rfixlen);
		mod= (i-(lm*quplen*sudakovlen*modellen*rfixlen)-(q*sudakovlen*modellen*rfixlen)-(s*modellen*rfixlen))/(rfixlen);
		rf= (i-(lm*quplen*sudakovlen*modellen*rfixlen)-(q*sudakovlen*modellen*rfixlen)-(s*modellen*rfixlen)-( mod*rfixlen));
		
		sprintf(command ,"mkdir %s/Mass%s-Qup%s-Model%s-Sud%s-rfix%s",dir,lmass[lm],qup[q],model[mod],sudakov[s],rfix[rf]);
		printf("%s\n",command );
		system(command);	
		
		sprintf(filenames[i],"%s/Mass%s-Qup%s-Model%s-Sud%s-rfix%s/control.h",dir,lmass[lm],qup[q],model[mod],sudakov[s],rfix[rf]);
		printf("%s\n",filenames[i]);	
		
		file=fopen(filenames[i],"w");
		if(file==NULL){
		 	printf("error opening file.");
		 	return 1;
		 }
		 if(strcmp(lmass[lm],"Default")!=0){
		 	fprintf(file,"#define MASS_L2 %s\n",lmass[lm]);
		 	fprintf(file,"#define MASS_S2 %s\n",lmass[lm]);
		 }
		 if(strcmp(qup[q],"Default")!=0){
		 	fprintf(file,"#define Q2_MAX %s.0\n",qup[q]);		 	
		 }
		 if(strcmp(model[mod],"Default")!=0){
		 	fprintf(file,"#define MODEL %s\n",model[mod]);		 	
		 }
		 if(strcmp(sudakov[s],"Default")!=0){
		 	fprintf(file,"#define SUDAKOV %s\n",sudakov[s]);		 	
		 }
		 if(strcmp(rfix[rf],"Default")!=0){
		 	fprintf(file,"#define R_FIX %s\n",rfix[rf]);		 	
		 }
		 fclose(file);
	}
	
	

	return 0;

}
