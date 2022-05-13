//#include<iostream>
//#include<fstream>
//#include<string>
//#include<vector>
//#include<cmath>
#include<stdio.h>
#include<math.h>


void simps5(double(*func)(double),double xmin,double xmax,double *res,double *err){
	//do simpsons approximation with n=2.
	// number of terms in the sum =2n+1.
	int n=2;

	double eval[section] ;
	double points[section];
	double step=(xmax-xmin)/(2*n);

	for(int i=0;i<(2*n);i++){
		points[i]=xmin+(i*step);
	}

	for(int i=0;i<(2*n);i++){
		eval[i]=(*func)(xmin+( points[i]) );
	}
	
	//evaluate uncertainty, see Bronstein et al .
	double diff_x[2*n];
	double diff_y[2*n];
	double fourth_d=0.0;//largest abs of fourth derivative...
	for(int i=0;i<4;i++){// i= number of deriv.
		for(int j=0;j<(2*n-i);j++){// then take gradient from sampled points.
			diff_x[j]=(points[j+1]+points[j])/2.0;//mid point;
			diff_y[j]=(eval[j+1]-eval[j])/(points[j+1]-points[j]);//grad
			
			if(i==3){//find max
				if( fourth_d < abs(diff_y[j]) ){
					fourth_d=abs(diff_y[j]);
				}
			}

		}
	}
	
	*err=fourth_d*(xmax-xmin)*pow(step,4)/180.0;

	*res=eva[0]+eval[section-1];


	for(i=1;i<=n ;i++){
		 *res+=4*eval[2*i-1];
	}
	for(i=1; i<n ;i++){
		*res+=2*eval[2*i];
	}
};


void simps_sector(double( *func)(double) ,double xmin, double xmax,double *res_arr,double *err_arr){
	int sector_number=sizeof(res_arr)/sizeof(*res_arr);

	double values[sector_number+1];//where the results for individual results go
	double errors[sector_number+1];//and correponding errors
	
	double step = (xmax-xmin)/sector_number;
	for(int i=0;i<sector_number; i++){
		
		simps5(&funcs, xmin+(step*i),xmin+(step*(i+1)), values+i, errors+i );
		fprintf(stdout, "Value: %f \t Error: %f \n", *(values+i),*(errors+i);
	}
	



}


