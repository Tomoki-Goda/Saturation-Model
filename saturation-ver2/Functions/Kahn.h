#ifndef KAHN
#define KAHN
//static const int N=3;
//struct Accumulator{
//	double *accum=(double*)malloc(N);
//} 

static int Kahn_init(double *accum,int N){
	for(int i=0;i<N;i++){
		accum[i]=0;
	}
	return 0; 
}
static double Kahn_total(double sum, double * accum,int N){
	double total=sum; 
	for(int i=0;i<N;i++){
		total+=accum[i];
	}
	//printf("total=%.5e\n",total);
	return total; 
}

static int plus(double *a,double *b){
	double c;
	double	sum,large, small;//=*a+*b;
	
	if(*a>*b){
		large=*a;
		small=*b;
	}else{
		large=*b;
		small=*a;
	}
	sum=large+small;
	c=sum-large;
	c=c-small;
	*a=sum;
	*b=c;
}
static int extract_exp(double x){
	int y;
	frexp(x,&y);
	return y;
}
static int accum_sort(double *accum,int len){
	int pos=0;
	while(1){
		if(pos+1==len){
			break;
		}
		//if(fabs(accum[pos])>fabs(accum[pos+1] )){
		//if(( (extract_exp(accum[pos])/5) > (extract_exp(accum[pos+1])/5) )|| accum[pos+1]==0.0 ){//5 is arbitrary, need two elements to be sufficiently different scales
		if(accum[pos]+accum[pos+1]==accum[pos]){
			pos++;
			continue;
		}else{
			plus(accum+pos,accum+pos+1);
			(pos==0)?(pos++):(pos--);
			continue;
		}
	}
	return 0;
}
static double Kahn(double a,double b,double *accum,int N){
	double accum_tmp=b;
	double sum=a;
	plus(&sum,&accum_tmp);
	int flag=1;
	while(flag){
		//printf("flag: %d\n",flag);
		for(int i=0;i<N;i++){
			plus(accum+i, &accum_tmp);
			//printf("accum: %.2e\t",accum_tmp);
			if(accum_tmp==0.0){
				flag=0;	
				break;
			}
			if(i+1==N&&flag==1){
				accum_sort(accum,N);
				i=0;
				flag=0;				
			}
		}
		//printf("\n");
	}
	return sum;
}

static double Kahn_list_sum(double* list, int len){
	double accum[3]={0};
	double sum=0;
	for(int i=0;i<len;i++){
		sum=Kahn(sum,list[i],accum,3);
	}
	sum=Kahn_total(sum,accum,3);
	return sum;
}



#endif
