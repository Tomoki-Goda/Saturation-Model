#ifndef KAHN_SUM_LOADED
#define KAHN_SUM_LOADED 1

#include<cmath>
#include<iostream>

///////////////////////////////  Usage  //////////////////////////////////////////
//  Create 
//  	Kahn sum
//  initialize
// 	sum.zero( a )
//  where a is a number of accumulator. 2 is usually enough...	
//  (name can be anything, but we are computing sum...)
//  use += operator to add,
//  	sum+= double ;
//  and total() to get total  	
//  	sum.total(); 
//////////////////////////////////////////////////////////////////////////////////


class Kahn{
	public:
		Kahn(int length){	
			len=length;
			accum=(double *)calloc(len,sizeof(double));
			//zero();
			sum=0;
		}
		//Kahn(void){	
		//	len=2;
		//	accum=(double *)malloc(len*sizeof(double));
		//}
		~Kahn(){
			free(accum);
		}
	private:
		int len;
		double *accum;
		double sum=0;
		
		int plus(double a, double b, double*sum, double *accum)const{
			double small=0,large=0;
			//double tot,acc;
			
			if(fabs(b)>fabs(a)){
				large=b;
				small=a;			
			}else{
				large=a;
				small=b;
			}
			(*sum)=large+small;
			*accum=small-( (*sum)-large );//be careful when using optimization 
			//(*accum)=small-(*accum);
			return 0;
		}
		int accum_sum(double val)const{
			double accum_tmp=val;
			for (int i=0;i<len;i++){
				
				plus(accum_tmp, accum[i], accum+i, &accum_tmp);
				//printf("%.5e \t ",accum[i]);
				if(accum_tmp==0.0){
					//printf("skip %d %e \n",i,accum_tmp );
					break;
				}
				if(i+1==len/*&& accum_tmp!=0.0*/){
					for(int j=0;j<len-1;j++){
						merge_accumulator(accum,len);
					}
				}
				
				
			}
			
			//printf("%.5e \t ",accum[len-1]);
			//printf("\n ");
			return 0;
		}
		
		double temp_accum=0;
		
		int merge_accumulator(double* accum ,int len)const{
			//double new_accum[len]={0};
			//Almost sorting algorithm.
			int pos=0;
			//for (int i=0;i<len;i++){
			//		printf("%.5e \t ",accum[i]);
			//	}
			//printf("\n ");
			while(1){
				//for (int i=0;i<len;i++){
				//		printf("%.5e \t ",accum[i]);
				//	}
				//printf("\n ");
				if(pos==len-1){
					break;
				}
				if(accum[pos]-accum[pos+1]==accum[pos]){
					//already good;
					pos+=1;
					continue;					
				}
				plus(accum[pos],accum[pos+1],accum+pos,accum+pos+1);
				if(pos==0){
					pos++;
				}else{
					pos-=1;
				}
				
			}	
			return 0;
		}
	public:	
		int zero(){
			sum=0;
			for(int i=0;i<len;i++){
				accum[i]=0;
			}
			return 0;
		}
		double total(int a){
			double val=total();
			zero();
			return(val);
		}

		double total(){
			double total=sum;
			double min=fabs(accum[0]);
			merge_accumulator(accum,len);
			for(int i=0;i<len;i++){
				if(min>fabs(accum[i])){
					min=fabs(accum[i]);
				}
				total+=accum[i];
			}	
			//if(min!=0.0){
			if(min>fabs(1.0e-15*total)){
				//merge_accumulator(accum,len);
				printf("Kahn:: Accumulator full: %.5e +- %.5e\n",total, min);
				for (int i=0;i<len;i++){
					printf("%.5e \t ",accum[i]);
				}
				printf("\n ");
			}
			return total;
		}
		
		double operator+=(double val){
			plus(sum,val,&sum,&temp_accum);
			accum_sum(temp_accum);
			return(sum);
		}			
	
};


//Test 
/*
int main(){
	double numbers[15]={
		-1.0e+28, -1.0e+18,1, 1.0e-16,1.0e+18,-1,1.0e+28,1.0e-30,
	 	-1.0e+18,-1.0e+28,1,-1,1.0e+28,-1.0e-16,1.0e+18};
	double sum=0;
	Kahn sum2(3);
	sum2.zero();
	
	for(int i=0;i<15;i++){
		sum2+=numbers[i];
		sum+=numbers[i];
		//printf("%.5e\n",sum2.total() );
	}

	printf("Result : Regular = %.5e\t Kahn = %.5e\n",sum,sum2.total() );
	
	return(0);	
}
*/
#endif
