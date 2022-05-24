#include<stdio.h>
#include<math.h>
#include"./dipole-cross-section.h"
#include<time.h>
//#include<file.h>

//FILE* outfile=fopen("./plots.m","w");

int main(){
	FILE* outfile=fopen("./plots.m","w");

	double par[7]={1.0,0.24,0.5e-4,1.26,4.0,0.2,0.8  };

	double Qs2=10;
	double x=1.0e-3;
	
	double r, rmin, rmax;
	rmax=10.0;
	rmin=0.001;
	unsigned point_n=1000;
	double pos, val;
	double step=( (double)(rmax-rmin) )/point_n;
	fprintf(outfile,"{");
	clock_t total;
	clock_t t;
	for(unsigned i=0;i<=point_n;i++){
		t=clock();
		pos=rmin+i*step;
		val=sigma(pos, x, Qs2, par ,2 , 'l');
		val+=sigma(pos, x, Qs2, par ,2 , 's');
		val+=sigma(pos, x, Qs2, par ,2 , 'c');
		val+=sigma(pos, x, Qs2, par ,2 , 'b');
		t-=clock();
		total-=t;

		fprintf(outfile, "{%f,%f},",pos, val);
		
	}
	printf("%f seconds\n",((float)total )/CLOCKS_PER_SEC);
	fprintf(outfile,"\b}");
	fclose(outfile);
	


}
