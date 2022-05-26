#include<stdio.h>
#include<math.h>

//extern "C" double simps_(double * , double *,
  //         double *,double *,double* ,double(*)(double*), double *  ,double *,double* ,double *);


void simpson1d(double(*function)(double ,double**), double ** par ,double min, double max,double *res){
	unsigned n=50;//2*n+1 terms in the sum
	double step= (max-min)/(2*n);
	double result=0;

	result+=(*function)(min,par)+(*function)(max,par);

	for(unsigned i=1;i <= n;i++){	
		(i!=n)?(  result+=4*(*function)(min+(2*i-1)*step ,par)+2*(*function)(min+2*i*step, par)  ) :( result+=4*(*function)(min+(2*i-1)*step ,par) );
	}
	result*=step /3;
	//printf("%.4e +/- N/A\n",result);
	*res=result;
	
}

void simpson1dA(double(*function)(double ,double**), double ** par ,double min, double max,unsigned n, double *res){
	//unsigned n=50;//2*n+1 terms in the sum
	double values[2*n+1];
	double diff[2*n];
	
	double step= (max-min)/(2*n);
	double result=0;
	double maxd4=0;

	for(unsigned i=0;i <= 2*n;i++){
		*(values+i)=(*function)(min+i*step,par);	
	}
	
	for(unsigned j=0;j<(2*n );j++){
			*(diff+j)=(*(values+j+1))-(*(values+j));
	}
	for(unsigned i=1;i<4;i++){
		for(unsigned j=0;j<(2*n -i);j++){
			*(diff+j)=(*(diff+j+1))-(*(diff+j));
			if((i==3 )&&( maxd4<fabs(*(diff+j))) ){ 
				maxd4=fabs(*(diff+j));
			}
		}
	
	}
	//max4d*=pow(step,-4);
	double error= maxd4*fabs(max-min)/180;
	for(unsigned i=0;i <= 2*n;i++){
		if((i==0)||(i==(2*n))){
		
		result+=(*(values+i));
		}else if((2*(i/2))==i){
		
		result+=2*(*(values+i));
		}else{
		
		result+=4*(*(values+i));
		}
	}
	result*=step /3;
	if(((error/result)>0.5)&& (error>1.0e-3)){
	printf("%.4e \t+/- %.4e\n",result,error);
	}
	*res=result;
	
}

	
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////	

unsigned simpcoeff(unsigned * pos,unsigned* len, unsigned dim){
	unsigned ctr=1;
	for(unsigned i=0;i<dim;i++){
		
		if((*(pos+i)==0)||(*(pos+i)==(*(len+i)) ) ){
		}else{
			((*(pos+i)/2)*2 ==*(pos+1))? (ctr*=2):(ctr*=4 ); 
		}

	}
	return( ctr);

}

void simpson2d(double(*function)(double,double, double**), double** par ,double xmin, double xmax,double ymin,double ymax,double *res){
        unsigned xn=50;//2*n+1 terms in the sum
	unsigned yn=50;

        double xstep= (xmax-xmin)/(2*xn);
	double ystep= (ymax-ymin)/(2*yn);

        double result=0;
        //double func=(*function);
	unsigned pos[2] ,len[2];
	(*len)=2*xn;
	*(len+1)=2*yn;
	

	unsigned coef;
	
        for(unsigned i=0;i<=2*xn;i++){
		for(unsigned j=0;j<=2*yn;j++){
			*pos=i;
			*(pos+1)=j;
			result+=( simpcoeff(pos,len,2) *(* function)(xmin+xstep*i,ymin+ystep*i, par));
					
		}
	}
	*res=(xstep*ystep/9)*result;
}




