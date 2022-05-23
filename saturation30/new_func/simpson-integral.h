#include<stdio.h>
#include<math.h>

//extern "C" double simps_(double * , double *,
  //         double *,double *,double* ,double(*)(double*), double *  ,double *,double* ,double *);

void simpson1d(double(*function)(double ,double**), double ** par ,double min, double max,double *res){
	if(fabs(max-min)<1.0e-15){
		*res=0.0;
		printf("min=max %.3e\n" ,*res);

	}else{

	unsigned n=25;//2*n+1 terms in the sum
	double step= (max-min)/(2*n);
	double result=0;
	//double func =(*function);

	result+=(*function)(min,par)+(*function)(max,par);
	//printf("init: %f\n",result);

	for(unsigned i=1;i <= n;i++){	
		(i!=n)?(  result+=4*(*function)(min+(2*i-1)*step ,par)+2*(*function)(min+2*i*step, par)  ) :( result+=4*(*function)(min+(2*i-1)*step ,par) );
		//printf("%f\n",result);
	}
	*res=result*step /3;
	//printf("%f %f %.3e\n" ,min,max,*res);
	}
}	

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
	*res=result;
}




