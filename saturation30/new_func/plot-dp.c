#include<stdio.h>
#include<math.h>
#include<time.h>

#include"./control.h"
#include"./constants.h"

#include"simpson-integral.h"
#include"dipole-cross-section.h"
#include"photon-wave-function.h"

#include"DIS-cross-section.h"

#include"cfortran.h"
#include"../minuit.h"
#include"./read-and-fit.h"
//#include<file.h>

//FILE* outfile=fopen("./plots.m","w");

int main(){
	FILE* testfile=fopen("./pwf.txt","w");
	double test_z,test_r,test_Q2, test_val;
	test_z=0.5;
	test_Q2=10.0;
	for(unsigned i =0;i<100;i++){
		test_r =pow(10,-3+((double)i)/100);
		test_val=psisq_f ( test_r,test_z,  test_Q2, 'l'); 
		//test_val+=psisq_f ( test_r,test_z,  test_Q2, 's');
		//test_val+=psisq_f ( test_r,test_z,  test_Q2, 'c');
		//test_val+=psisq_f ( test_r,test_z,  test_Q2, 'b');
		
		fprintf(testfile,"%f\t%f\n",test_r,test_val);
		
	}
	fclose(testfile);

	//FILE* outfile1=fopen("./dp-plots.txt","w");
	FILE* outfile2=fopen("./f2-plots.txt","w");

	double par[7]={23.0,0.28,3.0e-4,1.26,4.0,0.2,0.8  };

	double Qs2=10;
	double x=1.0e-3;
	
	double r, rmin, rmax;
	rmax=10.0;
	rmin=0.001;
	unsigned point_n=50;
	double pos, val1, val2;
	double step=( (double)(rmax-rmin) )/point_n;
	//fprintf(outfile1,"{");
	//fprintf(outfile2,"{");
	
	clock_t total;
	clock_t t;
	double Q2;
	for(unsigned i=0;i<=point_n;i++){
		t=clock();
		//pos=pow(10.0 , -3+((double)6*i)/point_n);
		
		Q2=pow(10.0 , -3+((double)6*i)/point_n);
		
		//val1=sigma_gbw(pos, mod_x(x,Qs2,'l'), Qs2, par );
		
		val2=sigma_DIS(x,Q2,1.0,par);
		
		
		t-=clock();
		total-=t;

		//fprintf(outfile1, "%f\t%f",pos, val1);
		//(i==point_n)?(fprintf(outfile1," ")):fprintf(outfile1,"\n") ;
		
		fprintf(outfile2, "%f\t%f",Q2, val2);
		(i==point_n)?(fprintf(outfile2," ")):fprintf(outfile2,"\n") ;
	}
	printf("%f seconds\n",((float)total )/CLOCKS_PER_SEC);
	
	//fclose(outfile1);
	fclose(outfile2);
	
//	system("math -i <\"plotdp.wls\"");
//	system("math -i <\"plotf2.wls\"");

}
