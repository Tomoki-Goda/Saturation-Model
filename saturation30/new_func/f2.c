#include<stdio.h>
#include<math.h>
#include<time.h>

#include"./control.h"
#include"./constants.h"


#include"simpson-integral.h"
#include"dipole-cross-section.h"
#include"photon-wave-function.h"

#include"DIS-cross-section.h"




int main(){
	fprintf(stdout,"started\n");

	FILE* file=fopen("./sigma.m", "w");
	//fprintf(file, "{");
	double x=1.0e-3;
		
	double par[7]={1.0,0.24,0.5e-4,1.26,4.0,0.2,0.8  };

	double res;
	double q2;
	

	fprintf(file, "{");
	clock_t time;
	for(unsigned i=0;i<=100;i++){
		q2=1.0e-3+10*i;
		time=clock();
		res=sigma_DIS(x,q2,1,par);
		time-=clock();
		fprintf(file, "{%f,%f}",q2,res);
		(i==100)?(fprintf(file,"}")):(fprintf(file,","));
		fprintf(stdout, "{%.3e,%.3e} : %.3e ,\n",q2,res, -((double)time)/CLOCKS_PER_SEC);

	}

//	fprintf(file,"\b}");
	fclose(file);
/*	
	file=fopen("./sigmaintegrand.m","w");
	fprintf(file,"{");
	double r,z,R;
	q2=100;
	z=0.5;	
	for(unsigned i=0;i<=1000;i++){
		//r=1.0e-3+(1.0e-2 *i);
		R=1.0e-4+1.0e-3 *i;
		r=R/(1-R);

		res=r *sigma_gbs(r,mod_x(x,q2,'l'),q2,par)* psisq_f(r,z,q2,'l');
		res*=1/pow(1-R,2);
		fprintf(file, "{%f,%f}",R,res);
		(i==1000)?(fprintf(file,"}")):(fprintf(file,","));

	}
	fclose(file);
*/
	return 0;
}

