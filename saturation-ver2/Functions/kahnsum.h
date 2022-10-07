#ifndef KAHN_SUM_LOADED
#define KAHN_SUM_LOADED 1
static void k_plus(double *sum, double add_term,double *counter){
	double lar,sma;
	if(fabs(*sum)>fabs(add_term)){
		lar=*sum;
		sma=add_term;
	}else{
		lar=add_term;
		sma=*sum;
	}
	*sum=lar+sma;
	double diff=*sum-lar;
	double coun=diff-sma;
	//return(sum);
}

/*static double k_group_sum(const double *arr, int len){
	double sum=0,count=0;
	for(int i=0;i<len;i++){
		k_plus(&sum,arr[i],&count);

	}
	return(sum+count);
}
*/

static double k_group_sum(const double *arr,int len){
	double sum50=0,sum25=0,sum0=0,sumM50=0,sumM25=0;
	double counter50=0, counter25=0,counter0=0,counterM25=0,counterM50=0;
	int ex;
	for(int i=0;i<len;i++){
		frexp(arr[i],&ex);
		if(ex>50){
			k_plus(&sum50,arr[i],&counter50);
		}else if(ex>25){
			k_plus(&sum25,arr[i],&counter25);
		}else if(ex>0){
			k_plus(&sum0,arr[i],&counter0);
		}else if(ex>-25){
			k_plus(&sumM25,arr[i],&counterM25);
		}else{
			k_plus(&sumM50,arr[i],&counterM50);
		}
	}
	//printf("%.5e,%.5e, %.5e, %.5e, %.5e\n ",counter50,counter25,counter0,counterM25,counterM50);
	sum50+=counter50;
	sum25+=counter25;
	sum0+=counter0;
	sumM25+=counterM25;
	sumM50+=counterM50;
	return(sum50+sum25+sum0+sumM25+sumM50);

}


//Test 
/*int main(){
	double numbers[5]={ -1.0e+18,1, 1.0e-18,-1,1.0e+18};
	double sum=0;
	double sum2=0;
	double counter=0;
	int exponent;
	double frac;
	for(int i=0;i<5;i++){
		frexp(numbers[i],&exponent);
		printf("%d\n",exponent);
		sum+=numbers[i];
		k_plus(&sum2,numbers[i],&counter);
	}

	sum2+=counter;
	double sum3=k_group_sum(numbers,5);
	printf("%.5e\t %.5e\t %.5e\n",sum,sum2,sum3);
	
	return(0);	
}*/
#endif
