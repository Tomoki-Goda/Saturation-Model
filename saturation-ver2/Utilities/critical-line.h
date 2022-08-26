
extern double SIGMA(double ,double, double ,double*,double *);
extern void approx_xg(double*);
extern int parameter(double*,double*,double*);
		

double sigma_func(double r, double** par){
	return(SIGMA(r,*(*(par)),*(*(par)+1),*(par+1),*(par+2))/(**(par+1))  );
}

double sigma_curve(double r, double** par){
	double x=**par;
	double Q2=*(*(par)+1);
	double *sigpar=*(par+1);
	double *sudpar= *(par+2);
	double dr=0.1*r;
	double curv=SIGMA(r+dr,x,Q2,sigpar,sudpar)+SIGMA(r-dr,x,Q2,sigpar,sudpar)- 2*SIGMA(r,x,Q2,sigpar,sudpar);
	return( curv/(dr*dr) );
}


double plot_curvature(double x,double Q2, double* sigpar, double* sudpar, FILE* file){
	double var[2];
	var[0]=x;
	var[1]=Q2;
	double * par[3];
	par[0]=var;
	par[1]=sigpar;
	par[2]=sudpar;
	double val,r;
	for(int i =0;i<1000;i++){
		r=pow(10,-1+3*((double)i)/1000);
		val=sigma_curve(r,par);
		fprintf(file,"%f\t%f\n",r,val);
	}
}

//double generate_point(double *par,double Q2,double x,char* filename){
double generate_points_crit(double x, double **param ){
	**param=x;

	double step=0.1;
	int max_n=100;
	double sample[max_n];
	double points[2]={0};
	double r;
	for(int i=0; i<=max_n;i++){
		r=pow(10,1-2*((double)i)/max_n);
		sample[i]=sigma_curve(r,param);
		if(sample[i]>1){
			for(int j =1 ;j<=i;j++){
				if((sample[i]*sample[i-j])<0){
					points[0]=pow(10,1-2*((double)i)/max_n);
					points[1]=pow(10,1-2*((double)(i-j))/max_n);
					break;
				}
			}
			break;
		}
	}
	if(points[0]==0||points[1]==0){
	printf("not found \n");
		getchar();
	}else{
	//	printf("do between %.3e %.3e\n",points[0],points[1]);
	}
///////////////////// find best point recursively /////////////////////////////////////////	
	double diffs[2];
//	double points[2];
	double aim=0.0; 

	double tolerance =1.0e-6;
	//points[0]=1.0e-1;
	//points[1]=20.0;
	double new_point;
	double new_diff;


	diffs[0]=sigma_curve(points[0],param);
	diffs[1]=sigma_curve(points[1],param);
	
	int error_flag=1;

	for(int i =0;i<50;i++){
		//new_point=(points[0]* fabs(diffs[1])+points[1]*fabs(diffs[0]))/( fabs(diffs[1])+fabs(diffs[0]) );
		new_point=(points[0] +points[1])/( 2 );
		new_diff=sigma_curve(new_point,param);
	//	if((i/10)*10==i){
	//		printf("%.3e :  %.5e %.5e %.5e\t %.5e %.5e %.5e\n",x,points[0], new_point, points[1],diffs[0], new_diff,diffs[1]);
	//	}
		if(fabs(new_diff)<tolerance){
			//printf("Found\n");
			error_flag=0;
		 	break;
		}else if((new_diff*diffs[0])<0 ){
			diffs[1]=new_diff;
			points[1]=new_point;
		}else if((new_diff*diffs[1])<0){
			diffs[0]=new_diff;
			points[0]=new_point;
		}else{
			printf("points: %.3e %.3e %.3e diff: %.3e %.3e %.3e\n", points[0],points[1],new_point,diffs[0],diffs[1],new_diff);
		}
	}
	if(error_flag!=0){
		printf("x: %.3e points: %.3e %.3e %.3e diff: %.3e %.3e %.3e\n",x, points[0],points[1],new_point,diffs[0],diffs[1],new_diff);
		getchar();
	}//else{
	//	printf("r = %f\n",new_point);
	//}
	return(2.0/(new_point*new_point));

}


