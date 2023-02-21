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
int N=256;
double wfull[129]={0}, whalf[65]={0}, x[129]={0};
std::string tag="unnamed";
int max_rec=7;
int InitDiv=1;
} CCIntegral;

//template <typename TYPE,typename args_type>static int fixed_cc(const CCIntegral & data,TYPE &func,args_type par,const double smin,const double smax,double&valfull,double &valhalf, double*arr){
template <typename TYPE,typename args_type>static int fixed_cc(const CCIntegral & data,TYPE &func,args_type par,const double smin,const double smax, Kahn &full, Kahn &half){
	const double *x16=data.x;
	const double *w16=data.wfull;
	const double *w8=data.whalf;
	const int N=data.N;
	double f[N+1];
	double arg;
	double scale=(smax-smin)/2;
	double mid=(smax+smin)/2;
	for(int i=0;i<N/2;i++){
		arg=scale*x16[i];
		f[2*i]=func(mid-arg,par);
		f[2*i+1]=func(mid+arg,par);
	}
	
	f[N]=func(mid,par);
	Kahn_clear(full);
	Kahn_clear(half);
	for(int i=0;i<N/2;i++){
		full+=f[2*i]*w16[i];
		full+=f[2*i+1]*w16[i];	
	}
	full+=f[N]*w16[N/2];
	for(int i=0;i<N/4;i++){
		half+=f[4*i]*w8[i];
		half+=f[4*i+1]*w8[i];
	}
	half+=f[N]*w8[N/4];
	full*=2*scale/N;
	half*=4*scale/N;
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
	double total=0;//,total2=0;
	const int N=data.N;
	double arr[N+1];
	Kahn accum=Kahn_init(3);
	Kahn accumfull=Kahn_init(3);
	Kahn accumhalf=Kahn_init(3);//, accum2=Kahn_init(3);
	double increase;
	smin=min;
	smax=(min-min/data.InitDiv)+max/data.InitDiv;//data.init_div;
	int licz=0,licztot=0 , counter=0;
			
	while(1){
		if(((max-min)-(smax-smin))==(max-min)||counter==MAX_RECURSION){
			smax=smin+2*scale;
			printf("Clenshaw_Curtis:: in \"%s\", evaluated %d times.\n",(data.tag).c_str(),counter );
			printf("sector size = %.3e\n [%.3e, %.3e] of [%.3e, %.3e] after %d / %d /%d\n",smax-smin,smin,smax,min,max, licz,licztot,MAX_RECURSION);
			
			//getchar();
			goto Error;
		}
		++counter;
		++licztot;
		scale=(smax-smin)/2;
		fixed_cc<TYPE,args_type>(data,func,par,smin,smax,accumfull,accumhalf);
		valfull=Kahn_total(accumfull);
		valhalf=Kahn_total(accumhalf);
#if DCLENSHAW_HH==1		
		if(not(std::isfinite(valfull)&&std::isfinite(valhalf))){
			printf("Clenshaw_Curtis:: in \"%s\" %.3e  %.3e encountered\n",(data.tag).c_str(),valfull,valhalf);
			goto Error;
		}
#endif
		accumhalf*=-1;
		accumhalf+=accumfull;
		//if(( fabs(valfull-valhalf)<eps*(fabs(valfull)) ) || (  fabs(valfull-valhalf)< fabs(smax-smin)*Aeps ) ){
		if(( fabs(Kahn_total(accumhalf))<eps*(fabs(valfull)) ) || (  fabs(Kahn_total(accumhalf))< fabs(smax-smin)*Aeps ) ){
			accum+=accumfull;
			++licz;
			counter=0;
			if(fabs(smax-max)==0.0){
				total=sign*Kahn_total(accum);
				Kahn_free(accum);
				Kahn_free(accumfull);
				Kahn_free(accumhalf);
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
		printf("valfull= %.3e , valhalf= %.3e  diff=%.3e\n",valfull,valhalf,Kahn_total(accumhalf));
		const double *x16=data.x;
		const double *w16=data.wfull;
		const double *w8=data.whalf;
		scale=(smax-smin)/2;
		double mid=(smax+smin)/2;
		for(int i=0;i<N/2;i++){
			arg=mid-scale*x16[i];
			printf("%.3e\t%.3e\n",arg,func(arg,par));
		}
		arg=mid;
		printf("%.3e\t%.3e\n",arg,func(arg,par));
		for(int i=0;i<N/2;i++){
			arg=mid+scale*x16[N/2-i-1];
			printf("%.3e\t%.3e\n",arg,func(arg,par));
		}
		printf("\n");
		Kahn_free(accum);
		Kahn_free(accumfull);
		Kahn_free(accumhalf);
		exit(1);
		//getchar();
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
	Kahn vec[N/4+1];
	double t[N/4+1];
	double c[N/2+1];
	if((N/8)*8!=N || N>256){
		printf("N=%d has to be multiple of 8, <= 128\n",N );
	}
	data.N=N;
	data.x[0]=1;
	data.x[1]=cos(PI/N);
	for(int i =2;i<N/2;i++){
		data.x[i]=cos(i*PI/N);
	}
	data.x[N/2]=0;
	
	//double *__restrict__ c=(double*)calloc((N/2+1),sizeof(double));
	c[0]=1;
	c[N/2]=1.0/(1-N*N);
	for(int i=1;i<N/2;i++){
		c[i]=2.0/(1-4*i*i);
	}
	//double *__restrict__ t=(double*)calloc((N/4+1),sizeof(double));
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
	//Kahn * vec=(Kahn*)malloc((N/4+1)*sizeof(Kahn));
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
			vec[j]*=t[j];
			accum+=vec[j];
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
			vec[2*j]*=t[2*j];
			accum+=vec[2*j];
			Kahn_clear(vec[2*j]);
		}
		data.whalf[i]=Kahn_total(accum);
		
		//data.whalf[i]+=accum;
		if(i==0){
			data.whalf[i]/=2;
		}
		//printf("%.10e\t",data.whalf[i]);
	}//printf("\n");
	
	//free(c);
	//free(t);
	Kahn_free(accum);
	for(int j=0;j<N/4+1;j++){
		Kahn_free(vec[j]);
	}
	//free(vec);
	//getchar();
	
	Kahn sum1=Kahn_init(2),sum2=Kahn_init(2);
	
	for(int i=0;i<N/2;i++){
		sum1+=2*data.wfull[i];
		if((i/2)*2==i){
			sum2+=2*data.whalf[i/2];
		}
	}
	sum1+=data.wfull[N/2];
	sum2+=data.whalf[N/4];
	double sum11=Kahn_total(sum1),sum22=Kahn_total(sum2);
	
	sum2*=-2;
	sum2+=sum1;
	
	if( fabs(Kahn_total(sum2)/(std::min(fabs(sum11),fabs(sum22))))>1.0e-15){
	//if( (sum1-2*sum2)!=0.0){
		printf("Factor array may be problematic, N = %d, %.3e - %.3e = %.3e\n",N,sum11,2*sum22,Kahn_total(sum2) );
		for(int i=0;i<N/2;i++){
			printf("%.3e ",data.wfull[i]);
			if((i/2)*2==i){
				printf("%.3e ",2*data.whalf[i/2]);
			}
			printf("\n");
		}
		getchar();
	}	
	Kahn_free(sum1);
	Kahn_free(sum2);
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
