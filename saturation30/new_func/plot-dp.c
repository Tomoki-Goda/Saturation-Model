#include<stdio.h>
#include<math.h>
#include"./dipole-cross-section.h"
//#include<file.h>

//FILE* outfile=fopen("./plots.m","w");

int main(){
	FILE* outfile=fopen("./plots.m","w");

	double par[5]={1.0,0.24,0.5e-4,1.26,5.0  };

	double Qs2=10;
	double x=1.0e-3;
	
	double r, rmin, rmax;
	rmax=10.0;
	rmin=0.001;
	unsigned point_n=1000;
	double pos, val;
	double step=( (double)(rmax-rmin) )/point_n;
	fprintf(outfile,"{");
	for(unsigned i=0;i<=point_n;i++){
		pos=rmin+i*step;
		val=sigma_gbs(pos, x, Qs2, par);

		fprintf(outfile, "{%f,%f},",pos, val);
		
	}
	fprintf(outfile,"\b}");
	fclose(outfile);
	


}
