#ifndef KAHN
#define KAHN 1 //change to one for testing
//static const int N=3;
//struct Accumulator{
//	double *accum=(double*)malloc(N);
//}
//#include<cassert>
typedef struct Kahn{double accum[5]; int N=0;} Kahn;
#ifndef KAHN_PREC 
	#define KAHN_PREC 1.0e-12
#endif

static void Kahn_clear(Kahn& kahn){
	for(int i=0;i<kahn.N;i++){
		kahn.accum[i]=0;
	}
}

static Kahn Kahn_init(int N){
	Kahn kahn;
	//kahn.accum=(double*)malloc(N*sizeof(double));
	kahn.N=N;
	Kahn_clear(kahn);
	return kahn;
}


static void Kahn_free(Kahn kahn){
	//free(kahn.accum);
}

static double Kahn_total(Kahn kahn){
	double total=0; 
	for(int i=0;i<kahn.N;i++){
		total+=kahn.accum[i];
	}
	return total; 
}
static void Kahn_show(Kahn kahn){
	printf(" {\t");
	for(int i=0;i<kahn.N;i++){
		printf("%.3e\t",kahn.accum[i]);
	}printf(" }\n");
}

static int plus(double &a,double &b){
	double c;
	double	sum,large, small;//=*a+*b;
	
	if(fabs(a)>fabs(b)){
		large=a;
		small=b;
	}else{
		large=b;
		small=a;
	}
	sum=large+small;
	c=sum-large;
	c=small-c;
	a=sum;
	b=c;
	return 0;
}

/*static int accum_sort(double (&accum)[5],int len){
	if(len==1){return 0;};
	int pos=0;
#if KAHN==1
	unsigned int counter=0;
#endif
	while(1){
#if KAHN==1
		if(++counter>10){printf("accum_sort:: Loop Error \n");
			printf("Ordered ?= {");
			for(int i=0;i<len-1;++i){ 
				printf("%d\t",accum[i]==(accum[i]+accum[i+1])); 
			}
			printf("}\t");
			
			printf("Accum = {");
			for(int i=0;i<len;++i){ 
				printf("%.3e\t",accum[i]); accum[i]=0;
			}
			printf("}\n");
			
			return 1;
		}
#endif
		if((accum[pos]+accum[pos+1])==accum[pos]){
			pos++;
			//continue;
		}else{
			//printf("sort %.2e %.2e %.3e\n", accum[pos],accum[pos+1], (accum[pos]+accum[pos+1])-accum[pos]);
			plus(accum[pos],accum[pos+1]);
			
			(pos==0)?(++pos):(--pos);
			continue;
		}
		if(pos+1==len){
			break;
		}
	}
	return 0;
}*/
static int orderedQ(double a, double b){
	//if a is much larger than b in magnitude, return 1
	int val;
	val=(fabs(b/a)<KAHN_PREC)?(1):(0);
	//val=(fabs(a+b)==fabs(a))?(1):(0);
	return val;
}
static int accum_sort(double (&accum)[5],int len){
	if(len==1){return 0;};
	int pos=0;
	for(int i=0;i<100;++i){
		if(orderedQ(accum[pos],accum[pos+1])!=1){
			if(i>3){printf("order : %.3e %.3e -> ",accum[pos],accum[pos+1]);}
			plus(accum[pos],accum[pos+1]);
			if(i>3){printf(" %.3e %.3e  \n",accum[pos],accum[pos+1]);}
			(pos==0)?(pos++):(pos--);			
		}else{
			++pos;
		}
		if(pos==len-1){
			break;
		}
		
	}
	return 0;
}
static int Kahn_Sum(Kahn& kahn, const double b){
	int flag=1;
	int N=kahn.N;
	double (&accum)[5]=kahn.accum;
	double accum_tmp=b;
	//accum[0]+=b;
	//return 0;
	
	unsigned int counter=0;
	while(flag){
		for(int i=0;i<N;i++){
			if(++counter>10){printf("Kahn_Sum::Loop Error \n");return 1;}
			plus(accum[i], accum_tmp);
			if(accum_tmp==0.0){
				flag=0;	
				break;
			}
			if(i==N-1){
				if(flag==1){
					accum_sort(accum,N);
					i=0;
					flag=0;	
				}			
			}
		}
	}
	return 0;
}
/*
static double Kahn_list_sum(double* list, int len){
	double accum[3]={0};
	double sum=0;
	for(int i=0;i<len;i++){
		sum=Kahn_Sum(sum,list[i],accum,3);
	}
	sum=Kahn_total(sum,accum,3);
	return sum;
}

*/

static Kahn& operator+=(Kahn& sum,const double a){
	Kahn_Sum(sum ,a);	
	return sum;
}
static Kahn& operator-=(Kahn& sum,const double a){
	Kahn_Sum(sum ,-a);	
	return sum;
}



static void Kahn_accum_sum(const Kahn& kahn1,Kahn& kahn2){
	if(kahn1.N!=kahn2.N){
		printf("incompatible accumulators. %d  %d \n ",kahn1.N,kahn2.N);
	}
	
	for(int i=0;i<kahn1.N;i++){
		kahn2+=kahn1.accum[i];
	}
}
static void Kahn_accum_diff(const Kahn& kahn1,Kahn& kahn2){
	if(kahn1.N!=kahn2.N){
		printf("incompatible accumulators. %d  %d \n ",kahn1.N,kahn2.N);
	}
	
	for(int i=0;i<kahn1.N;i++){
		kahn2-=kahn1.accum[i];
	}
}


static void Kahn_accum_times(const double a,Kahn& kahn2){
	for(int i=0;i<kahn2.N;i++){
		kahn2.accum[i]*=a;
	}
}

static Kahn& operator+=(Kahn& sum2,Kahn& sum){
	Kahn_accum_sum(sum ,sum2);	
	return sum2;
}
static Kahn& operator-=(Kahn& sum2,Kahn& sum){
	Kahn_accum_diff(sum ,sum2);	
	return sum2;
}

static Kahn& operator*=(Kahn& sum2,const double a){
	Kahn_accum_times(a,sum2);	
	return sum2;
}

static double Kahn_list_sum(double* list, int len){
	Kahn acc=Kahn_init(3);
	for(int i=0;i<len;i++){
		acc+=list[i];
	}
	return(Kahn_total(acc));
}
#endif
