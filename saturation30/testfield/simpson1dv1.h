////////////////////////////////////////////////////
//
// Tomoki Goda
// Simpson's approximation of 1 dimensional integration.
// May 2022
//
////////////////////////////////////////////////////
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

////////////////////////////////////////////////////
#define N 2 //simpsons approx with term 2*n+1. 
#define MAX_SECTOR 20
#define MIN_SECTOR 2 
////////////////////////////////////////////////////

void simpson_sum(double * values, double step, double* result ){
	double res=0;
	res=*(values)+*(values+(2*N));

	for(unsigned int i=1;i<=N;i++){
                res+=4* ( *( values+(2*i-1) ) );
        }

	for(unsigned int i=1;i< N;i++){
		res+=2* ( *( values+2*i ) );
	}
	res*=step/3;

	(*result)=res;
}
	
void simpson_error(double *points, double * eval, double * error ){
//evaluate uncertainty, see Bronstein et al .
	double diff_x[2*N];
	double diff_y[2*N];
	double fourth_d=0.0;//largest abs of fourth derivative...
	unsigned int len=2*N+1;
	for(unsigned int j=0;j<len-1 ;j++){// take gradient from sampled points.
                        *(diff_x + j)=(*(points+j+1)+*(points+j))/2;
                        *(diff_y + j)=(*(eval+ j+1)-*(eval+j))/(*(points+j+1)-*(points+j));
	}

	for(unsigned int i=1;i<4;i++){// i= number of deriv.
		for(unsigned int j=0;j<len-i-1 ;j++){// then take gradient from sampled points.
			*(diff_y + j)=(*(diff_y+ j+1)-*(diff_y+j))/(*(diff_x+j+1)-*(diff_x+j));
			*(diff_x + j)=(*(diff_x+j+1)+*(diff_x+j))/2;

			if(i==3){//find max
				if( fourth_d < fabs(*(diff_y+j)) ){
					fourth_d=fabs(*(diff_y+j));
				}
			}
		}
	}
	double range=(*(points+2*N)-*(points));
	double step=range/(2*N);
	//printf("step : %f \t range: %f\n", step, range);
	*error=fourth_d*(range )*pow(step,4)/180.0;
}

void simpson_evaluate(double (*function)(double ),double x_max,double x_min, double *x_val,double *y_val ){
	//this produces five points list for each sector.
	//x_val y_val are both pointer to the first-first element of 2d array for the sample points.
	unsigned int len =2*N*MIN_SECTOR+1;
	double eval[len];//array to store the evaluated points.
	double posi[len];//corresponding positions.
	
	double step=(x_max-x_min)/(len-1);

	for(unsigned int i=0;i < len;i++){
		*(posi+i)= x_min+i*step;
		*(eval+i)=(*function)(*(posi+i));
	}
	
	unsigned int sector_len=2*N+1;
	unsigned int position=0;
	for(unsigned int i=0; i<MIN_SECTOR;i++ ){
		for(unsigned int j=0;j<sector_len;j++){
			*(x_val+(i*sector_len)+j )=*(posi+position);
			*(y_val+(i*sector_len)+j )=*(eval+position);
			position++;
		}
		position--; 
	}


}

void sectors_summary(double *sect_res, double* sect_err,unsigned int sector_n, double *res, double *error, unsigned int * posit , double *max_err){
	//take results from sectors and add them and tell which sector has the largest uncertainty.
	unsigned int dpos=0;
	double max=0;
	(*res)=0;
	(*error)=0;
	for(unsigned int i=0;i<sector_n;i++){
		*res+=*(sect_res+i);
		*error+=*(sect_err +i);
		if( (*(sect_err+i)) > max ){
			dpos=i;
			max=*(sect_err+i);
		}
	}
	*posit=dpos;
	*max_err=max;
	//fprintf(stdout,"Current result: %.3e\t +- %.2e, ( %.1e percent ).\t%d sectors \n",*res,*error ,(100.0*(*error)/(*res)), sector_n );
	
}

void divide_sector(double(*func)(double),double* x_val,double* y_val, unsigned int position, unsigned int *sector_n, unsigned int * weight){
	//take list of 5 points for each sector. And at the give position, sector will be split,
	//resulting two sectors are stored in the original place and theother half will be appended to the list.
	//for this, resulting x_val and y_val are not in the right order after splitting. 
	//weight will be updated for the resultant sectors to account for the halving of the step size.
	unsigned int start=(2*N+1)*position;
	unsigned int appen=(2*N+1)*(*sector_n);

	double step=(*(x_val+start+2*N)- *(x_val+start ))/(2*2*N);
	
	for(unsigned int i=0; i<=N;i++){
	//first half is appended to the list;
		*(x_val+appen+2*i)=*(x_val+start + i);
		*(y_val+appen+2*i)=*(y_val+start + i);
		//new points.
		if(i!=N){
		*(x_val+appen+ 2*i+1)=(*(x_val+start + i))+step;	
		*(y_val+appen+ 2*i+1)=(*func)(*(x_val+appen+ 2*i+1) );
		}
	}	
	for(unsigned int i=N; i<=2*N;i++){
	//second half is stored in the original place 
		*(x_val+start+2*(i-N))=*(x_val+start + i);
		*(y_val+start+2*(i-N))=*(y_val+start + i);
		//new points.
		if(i!=2*N){
		*(x_val+start+ 2*(i-N)+1)=(*(x_val+start + i))+step;	
		*(y_val+start+ 2*(i-N)+1)=(*func)(*(x_val+start+ 2*(i-N)+1) );
		}
	}
	//now update weight since steps are 1/2 for divided sectors.
	//and sector_n in one larger.
	(*(weight+position))*=2;
	(*(weight+*sector_n))= (*(weight+position));
	(*sector_n)++;
			
}

void simpson_sector(double (*function)(double ),double x_max, double x_min,double aeps, double reps, double *res, double *error ){
	//main function
	//evaluate integral by dividing in small regions and use 2*N+1 term  simpson's approximation.
	//
	//double reps=0.001;
	//double aeps=0.001;
	
	*res=0.0;
	*error=0.0;
	unsigned int div_pos=0;
	double max_err=0.0;
	
	static unsigned int sector_len=2*N+1;
	double x_val[MAX_SECTOR][sector_len];
	double y_val[MAX_SECTOR][sector_len];
	double sect_res[MAX_SECTOR];
	double sect_err[MAX_SECTOR];
	double sect_rel[MAX_SECTOR];
	unsigned int sector_n=MIN_SECTOR;
	unsigned int weight[MAX_SECTOR];
	for(unsigned int i=0;i<MAX_SECTOR;i++){
		*(weight+i)=1;
	}
	///////////////initial evaluation ///////////////////////
	simpson_evaluate(function ,x_max, x_min,*x_val, *y_val );
	
	double step=(x_max-x_min)/(2*N*MIN_SECTOR);
	///////////////////////////////////////////////////////////////////////////////////
	//now we have x_val and y_val which are list of 5 points, 
	//then now convert them to the integral and error estimate.
	//and put them in sect_res and sect_err.
	//If quality is not good enough, add more sampling points by dividing sectors. 
	//until it reaches the max recursion or aimed error.
	for(unsigned j =0; j<=(MAX_SECTOR-MIN_SECTOR);j++){
		for(unsigned int i=0; i< sector_n;i++){
			simpson_error(((*x_val) +i*sector_len),((*y_val)+i*sector_len), (sect_err+i));
			simpson_sum(  ((*y_val) +i*sector_len),      step/(*(weight+i))        , (sect_res+i));
			*(sect_rel+i)=*(sect_err+i )/ (*(sect_res+i));
			//fprintf(stdout,"sector %d: %f\t +- %e \t Rel. %f\n",i+1, *(sect_res+i), *(sect_err+i),*(sect_rel+i) );		
			//printf(" weight: %d\n",(*(weight+i)) );
			//printf("\n");
		}
		//printf("\n");
		/////////////////////////////////////////////////////////////////////////////////
		
		sectors_summary(sect_res,sect_err, sector_n, res, error, &div_pos, &max_err);
		//fprintf(stdout,"Current result: %.4e\t +- %.2e, ( %.1e percent ).\t%d sectors \n",*res,*error ,(100.0*(*error)/(*res)), sector_n );
		//////////////////////////////////////////////////////////////////////////////////
		if((*error)<aeps){
			break;
		}
		if( ((*error)/(*res))<  reps  ){
			break;
		}
		//printf("Divide Sector %d\n", div_pos+1);
		////////////////sector division/////////////////////////	
		divide_sector(function,*x_val, *y_val, div_pos, &sector_n, weight);
		if(j>=(MAX_SECTOR-MIN_SECTOR)){
		        printf("reached limit.\n");
		}
	}
}


