double kahn_sum(const double *arr, int len){
//https://en.wikipedia.org/wiki/Kahan_summation_algorithm
	double sum=0 ;
	double c=0;
	
	double y,t;
	for(int i=0;i<len;i++){
		y=arr[i]-c;
		t=sum+y;
		c=(t-sum)-y;
		sum=t;
	}
	return sum;
}
double KBN_sum(const double *arr, int len){
//https://en.wikipedia.org/wiki/Kahan_summation_algorithm
	double sum=0 ;
	double c=0;
	
	double y,t;
	for(int i=0;i<len;i++){
		t=sum+arr[i];
		if(fabs(sum)>=fabs(arr[i])){
			c+=(sum-t)+arr[i];
		}else{
			c+=(arr[i]-t)+sum;
		}
		
		sum=t;
	}
	//printf("KBN accumulator: %.6e\n",c);
	return (sum+c);
}

