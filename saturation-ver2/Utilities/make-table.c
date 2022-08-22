#include<stdio.h>
#include<getopt.h>
//#include<unistd.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>


struct data{
	char sigma_0[20];
	
	char x_0[20];
	char lambda[20];
	char A_g[20];
	char lambda_g[20];
	char C1[20];
	char mu02[20];
	char C2[20];
	char mu022[20];
	char chisq[20];
	char dof[20];
	char chisqdof[20];
	

};

char testst[20]="-";

void read_write_data(struct data *data_struct, const char* name,char* str,const char* mode ){
	char* struc;
	char* dest;
	char* orig;
	

	if(strcmp(name,"sigma_0")==0){
		struc=data_struct->sigma_0;
	}else if(strcmp(name,"x_0")==0){
		struc=data_struct->x_0;
	}else if(strcmp(name,"lambda")==0){
		struc=data_struct->lambda;
	}else if(strcmp(name,"A_g")==0){
		struc=data_struct->A_g;
	}else if(strcmp(name,"lambda_g")==0){
		struc=data_struct->lambda_g;
	}else if(strcmp(name,"C")==0){
		struc=data_struct->C1;
	}else if(strcmp(name,"mu02")==0){
		struc=data_struct->mu02;
	}else if(strcmp(name,"C2")==0){
		struc=data_struct->C2;
	}else if(strcmp(name,"mu022")==0){
		struc=data_struct->mu022;
	}else if(strcmp(name,"chisq")==0){
		struc=data_struct->chisq;
	}else if(strcmp(name,"n_data-n_par")==0){
		struc=data_struct->dof;
	}else if(strcmp(name,"chisq/dof")==0){
		struc=data_struct->chisqdof;
	}else{
		printf("Unknown %s\n",name);
	}

	if(strcmp(mode,"w")==0){
		dest=struc;
		orig=str;
		//printf("write\n");
		//printf("%s\n",str);
	}else if(strcmp(mode,"r")==0){
		//printf("read\n");
		//printf("%s\n",struc);
		dest=str;
		orig=struc;
	}else{
		printf("Unknown mode. %s\n",mode);
	}
		strcpy(dest,orig);
	//printf("**** %s\n",struc);

}


//struct data name_arr[10];




int import_results(FILE * file,struct data *data_struct,struct data *error_struct, char vars[][20], int varlen ){
	char name[20],value[20],error[20];
	//fscanf(file, "%s\t%s",name,value);
	//printf("import\n");
	for(int j=0 ;j<100;j++){
		//printf("import\n");
		//fscanf(file, "%s\t%s\t%s",name,value,error);
		//printf("%s\n",name);	
		fscanf(file, "%s",name);
		//printf("name %s\n",name);	
		if(strcmp(name,"chisq")==0){
			//printf("Write\n");
			fscanf(file, "%s",value);
			read_write_data(data_struct,name,value,"w");
			fscanf(file, "%s",error);
			read_write_data(error_struct,name,error,"w");
			fscanf(file, "%s\t%s",name,value);	
			read_write_data(data_struct,name,value,"w");
			fscanf(file, "%s\t%s",name,value);	
			read_write_data(data_struct,name,value,"w");
			break;
		}
		for(int i=0;i<varlen;i++){
			if(strcmp(vars[i],name)==0){
				//printf("Write %s\n",name);
				fscanf(file, "%s",value);
				read_write_data(data_struct,name,value,"w");
				fscanf(file, "%s",error);
				read_write_data(error_struct,name,error,"w");
			}
		}
	}
}

void table_name(const char* key, char* name){
	if(strcmp(key,"sigma_0")==0){
		strcpy(name,"$\\sigma_0$");
	}else if(strcmp(key,"x_0")==0){
		strcpy(name,"$x_0 (10^{-4})$");
	}else if(strcmp(key,"lambda")==0){
		strcpy(name,"$\\lambda$");
	}else if(strcmp(key,"A_g")==0){
		strcpy(name,"$A_g$");
	}else if(strcmp(key,"lambda_g")==0){
		strcpy(name,"$\\lambda_g$");	
	}else if(strcmp(key,"C")==0){
		strcpy(name,"$C$");
	}else if(strcmp(key,"mu02")==0){
		strcpy(name,"$\\mu_0^2$");
	}else if(strcmp(key,"C2")==0){
		strcpy(name,"$C_S$");
	}else if(strcmp(key,"mu022")==0){
		strcpy(name,"$\\mu_{0S}^2$");
	}else if(strcmp(key,"chisq/dof")==0){
		strcpy(name,"$\\chi^2/dof$");
	}else if(strcmp(key,"chisq")==0){
		strcpy(name,"$\\chi^2$");
	}else{
		printf("Unknown\n");
	}
}

//char dir[500], mass[50], qup[50], model[50], sudakov[50] /*,rfix[50]*/,dum[100];
void interpret_dir(char* dir, char *mass, char * qup, char * model, char *sudakov){
		char* chptr;
		char* filenameptr[3];
		char dum[100];
		chptr=strtok(dir,"/");//split dir at /
		//printf("%s\n",chptr);
		while(chptr!=NULL){
			filenameptr[0]=filenameptr[1];
			filenameptr[1]=filenameptr[2];
			filenameptr[2]=chptr;
			chptr=strtok(NULL,"/");
		}
		sprintf(dir,"%s",filenameptr[0]);
		
		chptr=strtok(filenameptr[1],"-");
		sscanf(chptr,"%4s%s",dum,mass);
		chptr=strtok(NULL,"-");
		sscanf(chptr,"%3s%s",dum,qup);
		chptr=strtok(NULL,"-");
		sscanf(chptr,"%5s%s",dum,model);
		if((strcmp(model,"22")*strcmp(model,"2"))==0){
			strcpy(model,"$GBWS_{pert}$");
		}else if(strcmp(model,"0")==0){
			strcpy(model,"$GBW$");
		}else if(strcmp(model,"1")==0){
			strcpy(model,"$BGK$");
		}else if(strcmp(model,"3")==0){
			strcpy(model,"$BGKS_{pert}$");
		}else{
			printf("Unknown model %s\n",model);
		}		
		
		//chptr=strtok(NULL,"-");
		chptr=strtok(NULL,"/");
		sscanf(chptr,"%3s%s",dum,sudakov);
}

int main(int argc, char** argv){
	struct data data_arr[argc];
	struct data error_arr[argc];
	FILE* file;
	int sudakov_id, model_id;
	char shortopts[]="v:";
	struct option longopts[]={
		{"variables",required_argument,0,'v'},
		{"out",required_argument,0,'o'},
		{0,0,0,0}
	};
	int option_index=0;
	int signal=0;
	
	int directories[argc];
	int dir_ind=0;
	int output;

	char pars[argc][20];
	char *chptr=NULL;	
	int varlen=1;

	for(int i=1;i<argc;i++){

		signal=getopt_long_only(argc, argv,shortopts, longopts, &option_index);
		//printf("%d\n",signal);
		//if(signal==-1){
		//	break;
		//}

		switch(signal){
			case 'v':
				chptr=strtok(optarg," ");
				for(int j=0;j<10;j++){
					sprintf(pars[j],"%s",chptr);
					//printf("%s",chptr);
					chptr=strtok(NULL," ");	
					if(chptr==NULL){
						break;
					}
					varlen++;
				}
				i++;
				break;
			case 'o':
				i++;
				output=i;
				break;

			case -1:
				directories[dir_ind++]=i;
				//directories[++dir_ind]=0;
				//printf("%s\n",argv[i]);
				break;
			default:
				break;
		}
	}
	//printf("len= %d\tvarlen=%d\n",dir_ind,varlen);
	for(int i=0;i<dir_ind;i++){
		file=fopen(argv[directories[i]],"r");
		if(file==NULL){ printf("File not found:: %s\n",	argv[directories[i]]);}//else{printf("%dfile open\n",i);}
		import_results(file, &(data_arr[i]),&(error_arr[i]),pars,varlen);
	}
	//printf("out file \n");
	file=fopen(argv[output],"w");
	//printf("out file open\n");
	fprintf(file,"\\begin{table}\\resizebox{\\textwidth}{!}{\\begin{tabular}{|");
	
	for(int i=0;i<varlen+1;i++){
		fprintf(file,"c|");
	}
	fprintf(file,"}\n\\hline\n-");
	
	char parval[20];
	double parameter_value;
	for(int i=0;i<varlen;i++){
		table_name(pars[i],parval);
		fprintf(file,"&  %s ",parval);
		//printf(" $ %s $&",pars[i]);
	}
	
	struct data testdata= data_arr[0];
	//printf("chisq = %s\n", testdata.chisq); 
	

	fprintf(file,"\n\\\\\\hline\n");
	char qup[20], model[30], sudakov[20], mass[20];
	
	for(int i=0;i<dir_ind;i++){
		
		interpret_dir(argv[directories[i]], mass, qup,model, sudakov);
		fprintf(file,"%s, $m_l=%s$, $Q_{up}=%s$ ",model,mass,qup);
		for(int j=0;j<varlen;j++){
			//printf("%s ",pars[j]);
			read_write_data(&(data_arr[i]),pars[j],parval,"r");
			if(strcmp(parval,"")==0){
				strcpy(parval," - ");
			}
			parameter_value=atof(parval);
			//fprintf(file,"& %.3f ",atof(parval));
			if((fabs(parameter_value)<100)&&(fabs(parameter_value)>10)){
				fprintf(file,"& %.1f ",parameter_value);
			}
			else if((fabs(parameter_value)<10)&&(fabs(parameter_value)>1)){
				fprintf(file,"& %.2f ",parameter_value);
			}else if((fabs(parameter_value)<1)&&(fabs(parameter_value)>0.1)){
				fprintf(file,"& %.3f ",parameter_value);
			}else if((fabs(parameter_value)<0.1)&&(fabs(parameter_value)>0.01)){
				fprintf(file,"& %.4f ",parameter_value);
			}else{
				fprintf(file,"& %.3e ",parameter_value);
			}
		}
		fprintf(file,"\\\\ \\hline\n");
	}
	fprintf(file,"\\end{tabular} } \\end{table}\n");
	fclose(file);	

}




































