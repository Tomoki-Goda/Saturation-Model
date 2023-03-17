#include<iostream>
#include<getopt.h>
#include<string>
#include<fstream>
#include<stdlib.h>
#include<cmath>
#include<map>
#include<vector>
#include<stdio.h>
#include<cstring>
//#include <readline/readline.h>

typedef struct{
	std::string name;
	double value;
	double error;
} data;

int convert(char* str,double &val){
	val=atof(str);
	//printf("%s -> %f\n",str,val);
	return 0;
}

int convert(char* str, std::string& val){
	val=str;
	//printf("%s -> %s\n",str,val.c_str());
	return 0;
}

template<typename T >int store_data(char **lineptr ,T & slot){
	char *line=*lineptr;
	char letter=0, word[100];
	int j=0;
	while(1){
		letter=*(line);
		if(letter=='\0'){
			break;
		}
	//	printf("%c",letter);
		++line;
		if(letter==' '||letter=='\t'||letter=='\n' ){
			if(j==0){
				continue;
			}else{
				break;
			}
		}
		word[j++]=letter;
	}
//	printf("\n");
	word[j]='\0';
	convert(word,slot);
	*lineptr=line;
	return 0;
}
data store_data(char* line){
	data result={" ",0,0};
	store_data<std::string>(&line, result.name);
	store_data<double>(&line, result.value);
	store_data<double>(&line, result.error);
	
	return result;	
}
int read_file(std::fstream& file, std::map<std::string,double> &result, std::map<std::string,double> &error ){
	char c[100];
	size_t len;
	int pos=-1;
	data d;
	int count=0;
	while((!(file.eof()) )&&file ){
		if(++count>25){
			printf("we have problem\n");
			exit(1);
		}
		file.getline(c,100,'\n');
		d=store_data(c);	
		result.insert({d.name,d.value});
		error.insert({d.name,d.error});
		//printf("%s %.3e %.3e\n",d.name.c_str(),d.value,d.error);
	}
	return 0;
}


int split(std::string in, std::vector<std::string>& out){
	std::string::size_type pos1=0,pos2=0;
	std::string sub;
	while(1){
		while(in[pos1]==' '||in[pos1]=='\t'){
			++pos1;
		}
		pos2=in.find(",",pos1);
		if(std::string::npos==pos2){
			sub=in.substr(pos1,in.length()-pos1);
		}else{
			sub=in.substr(pos1,pos2-pos1);
		}
		
		std::cout<<pos1<<" "<<pos2<<"  "<<sub<<std::endl;
		//getchar();
		pos1=pos2+1;
		out.push_back(sub);
		if(pos2==std::string::npos){
			break;	
		}
		
	}
	return 0;
}

int file_exists(std::string name){
		//printf("checking %s \n",name.c_str());
		FILE* file=fopen(name.c_str(),"r");
		//int val=0;
		int val;//=(file==NULL)?0:1;
		if(file==NULL){
			val= 0;
		}else{
			val= 1;
			fclose(file);
		}
		return val;
}

/////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////

int main(int argc, char** argv){
	if((strcmp(argv[1],"--help")*strcmp(argv[1],"-h"))==0){
		printf("cmd \"par1-name, par2-name ...\" \"fit1-label, fit2-label, ... \" result1.txt result2.txt ... output-file.tex \n");
		printf(" column_names (sigma_0 etc.) row_names (GBW+sud etc. ) result_files tex_file\n"  );
		exit(0);
	}
	std::vector<std::string> names;
	std::vector<std::string> labels;
	std::map<std::string,std::string> print_names={
		{"sigma_0","$\\sigma_0 $ [mb]"},
		{"lambda","$\\lambda$"},
		{"x_0","$x_0 \\left(10^{-4}\\right)$"},
		{"A_g","$A_g$"},
		{"lambda_g","$\\lambda_g$"},
		{"mu02","$\\mu_0^2 \\left[\\mathrm{GeV^2}\\right]$"},
		{"C1","$C$"},
		{"mu102","$\\mu_0^2 \\left[\\mathrm{GeV^2}\\right]$"},
		{"C2","$C$"},
		{"mu202","$\\mu_0^2 \\left[\\mathrm{GeV^2}\\right]$"},
		{"chisq/dof","$\\chi^2/\\mathrm{dof}$"},
		{"Flag","Flag"}
      	};
	std::string name;
	std::string label;
	//char in[250];
	std::fstream file;
	FILE *	out;
	std::map<std::string,double> result[argc-3], error[argc-3];
	
	split(argv[1],names);
	split(argv[2],labels);

	char ans;
	
	if(file_exists(argv[argc-1])){
		printf("%s exists. Over-write? [y/n]\n",argv[argc-1]);
		scanf("%c",&ans);
		switch(ans){
			case 'y':
				break;
			case 'n':
				exit(0);
				break;
			default: 
				printf("unknown answer %c\n",ans);
				exit(1);	
			}
	}
	out=fopen(argv[argc-1],"w");
	fprintf(out,"\\begin{tabular}{|");
	for(int i=0;i<=names.size();++i){
		fprintf(out,"c|");
	}
	fprintf(out,"}\n\\hline");
	fprintf(out," - ");
	for(int i=0;i<names.size();++i){
		fprintf(out,"& %s ",print_names[names[i]].c_str());
	}
	fprintf(out,"\\\\\\hline \n");
	//FILE* filetest;
	for(int j=0;j<argc-4;++j){
		printf("%s\n",argv[3+j]);
		
		file.open(argv[3+j],std::fstream::in);
		read_file(file,result[j],error[j]);
		file.close();
		
		if(lrint(result[j]["Cov"])!=3||lrint(result[j]["Flag"])!=1){
			printf("Error Flags %f, %f, %f \n" ,result[j]["Flag"],result[j]["Error Flag"],result[j]["Cov"]);
			fprintf(out,"\\color{gray} (%1ld, %1ld) ",lrint(result[j]["Cov"]),lrint(result[j]["Flag"]));
		}
		fprintf(out,"{\\footnotesize %s} ", labels[j].c_str() );
		for(int i=0;i<names.size();++i){	
			name=names[i];	
			fprintf(out,"& %.3e",result[j][name]);
			printf("%.2e ",result[j][name]);
		}
		//fprintf(out,"\\\\\n");
		//for(int i=0;i<names.size();++i){	
		//	name=names[i];	
		//	if(error[j][name]!=0.0){
		//		printf(" +/- %.2e\t",error[j][name]);
		//		fprintf(out,"&\\hfill {\\tiny $\\pm$ %.3e}" ,error[j][name]);
		//	}else{
		//		fprintf(out,"& " );
		//	}
		//}
		fprintf(out,"\\\\\\hline \n");
		printf("\n");
	}
	fprintf(out,"\\end{tabular}");
	fclose(out);
	return 0;
}

