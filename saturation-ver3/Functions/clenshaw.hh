#ifndef DCLENSHAW_HH
#define DCLENSHAW_HH 1
#include<math.h>
#include<stdio.h>
#include<string>
#include"./Kahn.hh"
//#include<stdlib.h>

#ifndef PI
#define PI 3.141592653589793238462643383279502884197
#endif
typedef  struct{
int N=128;
double wfull[65]={0}, whalf[33]={0}, x[65]={0};
std::string tag="unnamed";
int max_rec=7;
int InitDiv=1;
} CCIntegral;

template <typename TYPE,typename args_type>static int fixed_cc(const CCIntegral & data,TYPE &func,args_type par,const double smin,const double smax,double&valfull,double &valhalf, double*arr){
	//const double (&x16)[]=data.x;
	//const double (&w16)[]=data.wfull;
	//const double (&w8)[]=data.whalf;
	const double *x16=data.x;
	const double *w16=data.wfull;
	const double *w8=data.whalf;
	const int N=data.N;
	//double *f=(double*)calloc(N+1,sizeof(double));
	double f[N+1];
	//double accum2[2]={0};
	Kahn val1=Kahn_init(3);
	//Kahn val2=Kahn_init(3);

	double arg;
	double scale=(smax-smin)/2;
	double mid=(smax+smin)/2;
	for(int i=0;i<N/2;i++){
		arg=scale*x16[i];
		f[2*i]=func(mid-arg,par);
		f[2*i+1]=func(mid+arg,par);
	}
	
	f[N]=func(mid,par);
	valfull=0;
	//Kahn_clear(val1);
	//Kahn_init(accum2,2);
	
	for(int i=0;i<N/2;i++){
		//val1+=f[2*i]*w16[i];
		//val1+=f[2*i+1]*w16[i];
		arr[2*i]=f[2*i]*w16[i];
		arr[2*i+1]=f[2*i+1]*w16[i];	
	}
	///val1+=f[N]*w16[N/2];
	arr[N]=f[N]*w16[N/2];
	Kahn_clear(val1);
	for(int i=0;i<N+1;i++){
		arr[i]*=2*scale/N;
		val1+=arr[i];	
	}
	valfull=Kahn_total(val1);
	valhalf=0;
	Kahn_clear(val1);
	
	//Kahn_init(accum2,2);
	for(int i=0;i<N/4;i++){
		val1+=f[4*i]*w8[i];
		val1+=f[4*i+1]*w8[i];
	}
	val1+=f[N]*w8[N/4];

	valhalf=Kahn_total(val1);
	//valfull*=2*scale/N;
	valhalf*=4*scale/N;
	Kahn_free(val1);
	//free(f);
	return 0;
}

template<typename TYPE,typename args_type>static double dclenshaw(const CCIntegral &data,TYPE &func, args_type par , const double a, const double b, const double eps, const double Aeps){
	if(a-b==0.0){
		return(0);
	}
	int MAX_RECURSION=data.max_rec;
//double dclenshaw(double(&func)(const double*, const void*),const void* args,const double a,const double b,const double eps){
	double sign, max,min;

	if(b>a){
		min=a;
		max=b;
		sign=1;
	}else if(a>b){
		max=a;
		min=b;
		sign=-1;
	}else{
		return(0);
	}
	double smin,smax;
	double valfull,valhalf;
	double scale;
	double arg;
	//double res;
	double total=0;//,total2=0;
	//double accum[2]={0},accum2[2]={0};
	const int N=data.N;
	double arr[N+1];
	//double *arr=(double*)calloc(N+1,sizeof(double));
	Kahn accum=Kahn_init(3);//, accum2=Kahn_init(3);
	//const int N=16;
	double increase;
	smin=min;
	//smax=max/data.div;
	smax=min+(max-min)/data.InitDiv;//data.init_div;
	int licz=0,licztot=0 , counter=0;
			
	while(1){
		if(((max-min)-(smax-smin))==(max-min)||counter==MAX_RECURSION){
			smax=smin+2*scale;
			printf("Clenshaw_Curtis:: in \"%s\", evaluated %d times.\n",(data.tag).c_str(),counter );
			printf("sector size = %.3e\n [%.3e, %.3e] of [%.3e, %.3e] after %d / %d /%d\n",smax-smin,smin,smax,min,max, licz,licztot,MAX_RECURSION);
			printf("valfull= %.3e , valhalf= %.3e  diff=%.3e\n",valfull,valhalf,valfull-valhalf);
			//getchar();
			goto Error;
		}
		++counter;
		++licztot;
		scale=(smax-smin)/2;
		fixed_cc<TYPE,args_type>(data,func,par,smin,smax,valfull,valhalf,arr);
#if DCLENSHAW_HH==1		
		if(not(std::isfinite(valfull)&&std::isfinite(valhalf))){
			printf("Clenshaw_Curtis:: in \"%s\" %.3e  %.3e encountered\n",(data.tag).c_str(),valfull,valhalf);
			goto Error;
		}
#endif
		//printf("[%.3e, %.3e] of [%.3e, %.3e] after %d / %d \n",smin,smax,min,max, licz,licztot);
		//printf("valfull= %.3e , valhalf= %.3e  diff=%.3e\n",valfull,valhalf,valfull-valhalf);
		
		if(( fabs(valfull-valhalf)<eps*(fabs(valfull)) ) || (  fabs(valfull-valhalf)<Aeps ) ){
			//PASS
			//accum+=valfull;
			//res=0;
			for(int i=0;i<N+1;i++){
				//res+=arr[i];
				accum+=arr[i];
			}
			//printf("%.3e %.3e %.3e %.3e \n",res,valfull,valhalf,res-valfull);
			//accum2+=valhalf;
			++licz;
			counter=0;
			if(fabs(smax-max)==0.0){
				//data.div/=licz;
				total=sign*Kahn_total(accum);
				//total2=sign*Kahn_total(accum2);
				Kahn_free(accum);
				//Kahn_free(accum2);

//				printf(" relative: %.3e / %.3e, %d  points /%d =%f \n", fabs(total2-total),total,data.N*licz, data.N*licztot,((double)licz)/licztot);
				//free(arr);
				return(total );
			}
			smin=smax;
			increase=(4*scale);
			smax=((max-(smin+increase)<(increase/2))?(max):(smin+increase));
		}else{
			//IMPROVE
			smax=smin+(scale/2);
		}

		if(((max-min)-(smax-smin))==(max-min)){
			printf("Clenshaw_Curtis:: in \"%s\", division exceeds limitation. in the domain [%.3e, %.3e] of [%.3e, %.3e] after %d / %d \n",(data.tag).c_str(), smin,smax,min,max, licz,licztot);						
			printf("valfull= %.3e , valhalf= %.3e \n",valfull,valhalf);
			//getchar();
			goto Error;
		}
		
	}
	
	Error:
		//const double (&x16)[]=data.x;
		//const double (&w16)[]=data.wfull;
		//const double (&w8)[]=data.whalf;
		//getchar();
		const double *x16=data.x;
		const double *w16=data.wfull;
		const double *w8=data.whalf;
		scale=(smax-smin)/2;
		double mid=(smax+smin)/2;
		//const int N=data.N;
		//getchar();
		for(int i=0;i<N/2;i++){
			arg=mid-scale*x16[i];
			//printf("f(%.3e) = %.3e\n",arg,func(arg,par));
			printf("%.3e\t%.3e\n",arg,func(arg,par));
			/*printf("%.3e\t%.3e\t%.3e\t%.3e",arg,func(arg,par), arr[2*i]*N/(2*scale*w16[i])  ,w16[i]);
			((i/2)*2==i)?
			(printf("\t%.3e\n",2*w8[i/2])):
			(printf("\n")) ;*/
		}
		arg=mid;
		//printf("f(%.3e) = %.3e\n",arg,func(arg,par));
		printf("%.3e\t%.3e\n",arg,func(arg,par));
		//printf("%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",arg,func(arg,par),arr[N]*N/(2*scale*w16[N/2])  ,w16[N/2],2*w8[N/4]);
		for(int i=0;i<N/2;i++){
			arg=mid+scale*x16[N/2-i-1];
			//printf("f(%.3e) = %.3e\n",arg,func(arg,par));
			printf("%.3e\t%.3e\n",arg,func(arg,par));/*
			printf("%.3e\t%.3e\t%.3e\t%.3e",arg,func(arg,par),arr[N-2*i-1]*N/(2*scale*w16[N/2-i-1]) ,w16[N/2-i-1]);
			((i/2)*2!=i)?
			(printf("\t%.3e\n",2*w8[N/4-i/2-1])):
			(printf("\n")) ;*/
			//printf("%.3e\t%.3e\t%.3e\n",arg, func(arg,par),w16[N/2-i-1]);
		}
		printf("\n");
		
		getchar();
		//free(arr);
		return 0;
	
}


inline int t_pos(int i,int j, int n){
	return(n/4-abs(n/2-i*j+((i*j)/n)*n));
}
inline int sign(int i){
	return( (i==0)?(1): (i/abs(i)));
}

CCIntegral CCprepare(const int N){
	CCIntegral data;
	
	if((N/8)*8!=N || N>128){
		printf("N=%d has to be multiple of 8, <= 128\n",N );
	}
	data.N=N;
	data.x[0]=1;
	data.x[1]=cos(PI/N);
	for(int i =2;i<N/2;i++){
		data.x[i]=cos(i*PI/N);
	}
	data.x[N/2]=0;
	
	double *__restrict__ c=(double*)calloc((N/2+1),sizeof(double));
	c[0]=1;
	c[N/2]=1.0/(1-N*N);
	for(int i=1;i<N/2;i++){
		c[i]=2.0/(1-4*i*i);
	}
	double *__restrict__ t=(double*)calloc((N/4+1),sizeof(double));
	//double accum[3]={0};
	Kahn accum=Kahn_init(4);
	
	for(int j=0;j<N/4+1;j++){
		if(j==0){
			t[j]=0;
		}else if(j==N/4){
			t[j]=-1;
		}else{
			t[j]=-sin(2*j*PI/N);
		}
		//printf("%.3e\t",t[j]);
		//t[i*(N/2+1)+j]=((i*j==0)?(1):(cos(2*i*j*PI/N)) );
	}//printf("\n");
	int pos;
	Kahn * vec=(Kahn*)malloc((N/4+1)*sizeof(Kahn));
	for(int j=0;j<N/4+1;j++){
		vec[j]=Kahn_init(4);
		Kahn_clear(vec[j]);
	}
	
	for(int i=0;i<N/2+1;i++){
		Kahn_clear(accum);
		data.wfull[i]=0;
		for(int j=0;j<N/2+1;j++){
			pos=t_pos(i,j,N);
			//accum+= sign(pos) * t[abs(pos)] *c[j];
			vec[abs(pos)]+=sign(pos)*c[j];
		//	printf("%.2e\t",t[abs(pos)]);
		}//printf("\n");
		for(int j=0;j<N/4+1;j++){
			Kahn_accum_times(t[j],vec[j]);
			Kahn_accum_sum(vec[j],accum);
			Kahn_clear(vec[j]);
		}
		
		data.wfull[i]=Kahn_total(accum);
		if(i==0){
			data.wfull[i]/=2;
		}
		//printf("%.10e\t",data.wfull[i]);
	}//printf("\n");
	
	//t[0]=1;
	for(int i=0;i<N/4+1;i++){
		Kahn_clear(accum);
		data.whalf[i]=0;
		for(int j=0;j<N/4+1;j++){
			//t[j]=cos(i*j*4*PI/N);
			pos=t_pos(i,j,N/2);
			vec[2*abs(pos)]+=sign(pos)*c[j];
			//accum+=sign(pos)*t[2*abs(pos)]*c[j];
		}
		for(int j=0;j<N/8+1;j++){
			Kahn_accum_times(t[2*j],vec[2*j]);
			Kahn_accum_sum(vec[2*j],accum);
			Kahn_clear(vec[2*j]);
		}
		data.whalf[i]=Kahn_total(accum);
		
		//data.whalf[i]+=accum;
		if(i==0){
			data.whalf[i]/=2;
		}
		//printf("%.10e\t",data.whalf[i]);
	}//printf("\n");
	
	free(c);
	free(t);
	Kahn_free(accum);
	for(int j=0;j<N/4+1;j++){
		Kahn_free(vec[j]);
	}
	free(vec);
	//getchar();
	double sum1=0, sum2=0;
	for(int i=0;i<N/2;i++){
		sum1+=2*data.wfull[i];
		if((i/2)*2==i){
			sum2+=2*data.whalf[i/2];
		}
	}
	sum1+=data.wfull[N/2];
	sum2+=data.whalf[N/4];
	if( fabs((sum1-2*sum2)/(std::min(fabs(sum1),fabs(sum2))))>1.0e-15){
	//if( (sum1-2*sum2)!=0.0){
		printf("Factor array may be problematic, N = %d : %.3e - %.3e = %.3e\n",N,sum1,2*sum2,sum1-2*sum2);
		for(int i=0;i<N/2;i++){
			printf("%.3e ",data.wfull[i]);
			if((i/2)*2==i){
				printf("%.3e ",2*data.whalf[i/2]);
			}
			printf("\n");
		}
		getchar();
	}	
	//getchar();
	return data;
}
CCIntegral CCprepare(const int N,const std::string &tag){
	CCIntegral data=CCprepare(N);
	data.tag=tag;
	
	return data;
}
CCIntegral CCprepare(const int N,const std::string &tag,int d){
	CCIntegral data=CCprepare(N);
	data.tag=tag;
	data.InitDiv=d;
	
	return data;
}
CCIntegral CCprepare(const int N,const std::string &tag,int d,int max){
	CCIntegral data=CCprepare(N);
	data.tag=tag;
	data.InitDiv=d;
	data.max_rec=max;
	
	return data;
}


#endif
