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
	
	
	char name[500];
	float value;
	float error;
	float chisq;
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
			printf("makelist error no file: %s\n",argv[i]);
			return 1;
		}
		fscanf(resfile,"%s\t%f",name, &value);//Qup
		//fprintf(outfile,"%.0f",value);
		char dir[500], mass[50], qup[50], model[50], sudakov[50] ,rfix[50],dum[100];
		char* chptr;
		char* filenameptr[3];

		chptr=strtok(argv[i],"/");//split dir at /
		//printf("%s\n",chptr);
		while(chptr!=NULL){
			filenameptr[0]=filenameptr[1];
			filenameptr[1]=filenameptr[2];
			filenameptr[2]=chptr;
			chptr=strtok(NULL,"/");
		}
		//files are usually model-name/Mass0.0-Qup650-.../result.txt 
		//printf("%s\n",filenameptr[0]);
		//printf("%s\n",filenameptr[1]);
		sprintf(dir,"%s",filenameptr[0]);
		
		chptr=strtok(filenameptr[1],"-");
		sscanf(chptr,"%4s%s",dum,mass);
		
		chptr=strtok(NULL,"-");
		sscanf(chptr,"%3s%s",dum,qup);
		chptr=strtok(NULL,"-");
		sscanf(chptr,"%5s%s",dum,model);
		chptr=strtok(NULL,"-");
		sscanf(chptr,"%3s%s",dum,sudakov);	
		chptr=strtok(NULL,"/");
		sscanf(chptr,"%4s%s",dum,rfix);
		
		if(i==1){
			if((strcmp(model,"0")*strcmp(model,"2")*strcmp(model,"22"))==0 ){
				fprintf(outfile,"\\begin{sidewaystable}\n");
				fprintf(outfile,"\\rowcolors{1}{}{lightgray}\n");
				fprintf(outfile,"\\resizebox{\\textwidth}{!}{\n");
				fprintf(outfile,"\\begin{tabular}{|c||c|c|c|c|c|c|c||c|c|}\n\\hline \n");
				fprintf(outfile,"-& $\\sigma_0$ & $\\lambda$ & $x_0 (10^{-4})$ & $C$ & $\\mu_0^2$ & $g_1$ & $g_2$ & $\\chi^2$ & Flag \\\\ \\hline \\hline \n");
			}
			if((strcmp(model,"1")*strcmp(model,"3"))==0 ){
				fprintf(outfile,"\\begin{sidewaystable}\n");
				fprintf(outfile,"\\rowcolors{1}{}{lightgray}\n");
				fprintf(outfile,"\\resizebox{\\textwidth}{!}{\n");
				fprintf(outfile,"\\begin{tabular}{|c||c|c|c|c|c|c|c|c|c||c|c|} \n\\hline \n");
				fprintf(outfile,"-& $\\sigma_0$ & $A_g$ & $\\lambda_g $ & $C$ & $\\mu_0^2$ & $C_s$ & $\\mu_{0s}^2$ & $g_1$ & $g_2$ & $\\chi^2$ & Flag \\\\ \\hline \\hline \n");
			}
		}
				
		//sscanf(argv[i],"./%s/Mass%s-Qup%s-Model%s-Sud%s-rfix%s/result.txt",dir,mass,qup,model,sudakov,rfix);
		//printf("./%s/Mass%s-Qup%s-Model%s-Sud%s-rfix%s/result.txt",dir,mass,qup,model,sudakov,rfix);
		//getchar();
		fprintf(outfile,"{\\footnotesize "); 
		//fprintf(outfile,"{\\tiny "); 
		if(strcmp(model,"0")==0){
			fprintf(outfile,"%s: GBW $m_l=%s$ $Q_{up}=%s$",dir,mass,qup);
		}else if(strcmp(model,"1")==0){
			fprintf(outfile,"%s: BGK $m_l=%s$ $Q_{up}=%s$",dir, mass,qup);
		}else if(strcmp(model,"2")==0||strcmp(model,"22")==0){
					if(strcmp(sudakov,"0")==0){
						fprintf(outfile,"%s: GBW(S) $m_l=%s$ $Q_{up}=%s$",dir, mass,qup);
					}else if(strcmp(sudakov,"1")==0){
						fprintf(outfile,"%s: $\\mathrm{GBWS_{pert}}$ $m_l=%s$ $Q_{up}=%s$",dir, mass,qup);
					}else if(strcmp(sudakov,"2")==0){
						fprintf(outfile,"%s: $\\mathrm{GBWS_{np}}$ $m_l=%s$ $Q_{up}=%s$",dir, mass,qup);
					}else{
						printf("invalid sudakov flag: %s\n",sudakov);
					}
					if(strcmp(rfix,"1")==0){
						fprintf(outfile," rfix");
					}
		}else if(strcmp(model,"3")==0){
					if(strcmp(sudakov,"0")==0){
						fprintf(outfile,"%s: BGK(S) $m_l=%s$ $Q_{up}=%s$",dir, mass,qup);
					}else if(strcmp(sudakov,"1")==0){
						fprintf(outfile,"%s: $\\mathrm{BGKS_{pert}}$ $m_l=%s$ $Q_{up}=%s$",dir, mass,qup);
					}else if(strcmp(sudakov,"2")==0){
						fprintf(outfile,"%s: $\\mathrm{BGKS_{np}}$ $m_l=%s$ $Q_{up}=%s$",dir, mass,qup);
					}else{
						printf("invalid sudakov flag: %s\n",sudakov);
					}		
					if(strcmp(rfix,"1")==0){
						fprintf(outfile," rfix");
					}
		}else{
			printf(" invalid model id : %s\n",model);
		}
		fprintf(outfile,"} ");
		
		int skip=0;
		if(strcmp(sudakov,"1")*strcmp(sudakov,"2")==0){
			skip=2;
		}

		//for(unsigned line =0; (line<(line_n)) ; line++){
		for(unsigned line =0; line<10 /* > max possible no. of parameter*/; line++){
			////////////////////////////////
			
			if(( strcmp(model,"0")*(abs( strcmp(model,"22")*strcmp(model,"2"))+abs(strcmp(sudakov,"0")) ))==0  ){
				if((3<=line)&&(line<=6)){
					fprintf(outfile," & %c ",'-');
					continue;
				}
			}
			if(( strcmp(model,"1")*( abs(strcmp(model,"3"))+abs(strcmp(sudakov,"0")) ) )==0  ){
				if((5<=line)&&(line<=8)){
					fprintf(outfile," & %c ",'-');
					continue;
				}
			}
			if(( abs(strcmp(model,"22")*strcmp(model,"2"))+abs(strcmp(sudakov,"1")))==0  ){
				if((5<=line)&&(line<=6)){
					fprintf(outfile," & %c ",'-');
					continue;
				}
			}
			if(( ( abs(strcmp(model,"3"))+abs(strcmp(sudakov,"1"))))==0  ){
				if((7<=line)&&(line<=8)){
					fprintf(outfile," & %c ",'-');
					continue;
				}
			}
			
			///////////////////////////////
			fscanf(resfile,"%s\t%f\t%f\n",name,&value, &error );
			//fprintf(outfile," & %.2e {\\tiny $\\pm$ %.2e } ",value, error);
			
			if(strcmp(name,"chisq")==0){
				printf("%d parameters \n",count);
				break;	
			}
			count++;
			fprintf(outfile," & %.2e {\\tiny $\\pm$ %.2e } ",value, error);
		}
		fprintf(outfile," &{\\footnotesize %.2e} {\\tiny $\\pm$ %.2e } ",value, error);
		//fscanf(resfile,"%s\t%f\t%f\n",name,&chisq, &error );
		//fprintf(outfile,"& %.2e  ",chisq);
		//fprintf(outfile,"& %.2e  ",value);
		double chisqerr=error;

		fscanf(resfile,"%s\t%d\n",name,&ndata);
		fprintf(outfile,"{\\footnotesize / %d }=",(ndata+skip));
		fprintf(outfile," %.2f  ",value/(ndata+skip));
			
		fscanf(resfile,"%s\t%f\n",name,&value);
		printf(	"%s\t%f\n",name,value);

		fscanf(resfile,"%s %s\t%d\n",name,name,&ndata);//error_flag
		fprintf(outfile,"&{\\footnotesize");
		if((ndata!=3)||( chisqerr>1.0e-4) ){
			fprintf(outfile," {\\color{red}%d} ",ndata);
		}else{
			fprintf(outfile," %d  ",ndata);
		}
		
		fprintf(outfile,"}\\\\ \n\\hline \n");


		fclose(resfile);	
	}
	fprintf(outfile,"\\end{tabular}\n");
	fprintf(outfile,"}\n");
	fprintf(outfile,"\\end{sidewaystable}\n");


	fclose(outfile);
	return 0;
}
