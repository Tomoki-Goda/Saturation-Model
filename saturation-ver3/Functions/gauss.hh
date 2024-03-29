#ifndef DGAUSS_H
#define DGAUSS_H
#include<math.h>
#include<stdio.h>
#include"./Kahn.hh"

static const double x10[]={
0.148874338981631211,	
0.433395394129247191,
0.679409568299024406,
0.865063366688984511,
0.973906528517171720	
};
static const double w10[]={
0.295524224714752870,
0.269266719309996355,
0.219086362515982044,
0.149451349150580593,
0.066671344308688138
};


static const double x20[]={
0.076526521133497333755,
0.227785851141645078080,
0.373706088715419560673,
0.510867001950827098004,
0.636053680726515025453,
0.746331906460150792614,
0.839116971822218823395,
0.912234428251325905868,
0.963971927277913791268,
0.993128599185094924786};

static double w20[]={
0.152753387130725850698,
0.149172986472603746788, 
0.142096109318382051329,
0.131688638449176626898,
0.118194531961518417312,	
0.101930119817240435037,
0.083276741576704748725,
0.062672048334109063570,
0.040601429800386941331,
0.017614007139152118312};


static const double x40[]={
0.038772417506050821933,	
0.116084070675255208483,	
0.192697580701371099716,	
0.268152185007253681141,	
0.341994090825758473007,	
0.413779204371605001525,	
0.483075801686178712909,	
0.549467125095128202076,	
0.612553889667980237953,	
0.671956684614179548379,	
0.727318255189927103281,	
0.778305651426519387695,	
0.824612230833311663196,	
0.865959503212259503821,	
0.902098806968874296728,	
0.932812808278676533361,	
0.957916819213791655805,	
0.977259949983774262663,	
0.990726238699457006453,	
0.998237709710559200350,	
};

static const double w40[]={
0.077505947978424811264,
0.077039818164247965588,
0.076110361900626242372,
0.074723169057968264200,
0.072886582395804059061,
0.070611647391286779696,
0.067912045815233903826,
0.064804013456601038075,
0.061306242492928939167,
0.057439769099391551367,
0.053227846983936824355,
0.048695807635072232061,
0.043870908185673271992,
0.038782167974472017640,
0.033460195282547847393,
0.027937006980023401099,
0.022245849194166957262,
0.016421058381907888713,
0.010498284531152813615,
0.004521277098533191258,
};
template<typename TYPE,typename args_type>static double dgauss40(TYPE &func, args_type param , const double smin, const double smax,Kahn& val40,Kahn& val20){
	//double accum[3]={0};
	//Kahn accum=Kahn_init(3);
	double x; 
	double scale=(smax-smin)/2;
	double mid=(smax+smin)/2;

	Kahn_clear(val20);
	Kahn_clear(val40);
	for(int i=0;i<10;i++){
		x=scale*x20[i];
		val20+=w20[i]*func(mid-x,param);
		val20+=w20[i]*func(mid+x,param);
	}
	//val20=Kahn_total(accum);
	//Kahn_clear(accum);
	//val40=0;
	for(int i=0;i<20;i++){
		x=scale*x40[i];
		val40+=w40[i]*func(mid-x,param);
		val40+=w40[i]*func(mid+x,param);
	}
	//val40=Kahn_total(accum);
	val20*=scale;
	val40*=scale;
	//Kahn_free(accum);

	return(Kahn_total(val40));
}
template<typename TYPE,typename args_type>static double dgauss20(TYPE &func, args_type param , const double smin, const double smax,Kahn& val20,Kahn& val10){
	//double accum[3]={0};
	//Kahn accum=Kahn_init(3);
	double x; 
	double scale=(smax-smin)/2;
	double mid=(smax+smin)/2;

	Kahn_clear(val10);
	Kahn_clear(val20);
	for(int i=0;i<5;i++){
		x=scale*x10[i];
		val10+=w10[i]*func(mid-x,param);
		val10+=w10[i]*func(mid+x,param);
	}
	//val20=Kahn_total(accum);
	//Kahn_clear(accum);
	//val40=0;
	for(int i=0;i<10;i++){
		x=scale*x20[i];
		val20+=w20[i]*func(mid-x,param);
		val20+=w20[i]*func(mid+x,param);
	}
	//val40=Kahn_total(accum);
	val10*=scale;
	val20*=scale;
	//Kahn_free(accum);

	return(Kahn_total(val20));
}



template<typename TYPE,typename args_type>static double dgauss(TYPE &func, args_type param , const double min, const double max, const double eps, const double epsabs){
	double smin,smax; //section
	int counter=0;
	int MAX_RECURSION=4;
	double val10=0, val20=0;
	double scale;
	double prec=1.0e-14;
	//double accum_total[3]={0};
	Kahn accum=Kahn_init(1);
	Kahn accum20=Kahn_init(1),accum10=Kahn_init(1);
	double error=0;
	double increase;
	double error_ratio;
	smax=max;
	smin=min;
	int licz=0,licztot=0;
	//double f[N/2+1];
	double f[11];
	if(min==max){
		return 0;
	}		
	while(1){
		if(counter==MAX_RECURSION){
			printf("dgauss::MAX_RECURSION\n");
			printf("[%.3e, %.3e] of [%.3e, %.3e] after %d / %d \n",smin,smax,min,max, licz,licztot);
			printf("valfull= %.3e , valhalf= %.3e  diff=%.3e\n",val20,val10,val20-val10);
			Kahn_free(accum);
			Kahn_free(accum10);
			Kahn_free(accum20);
			goto Error;
		}
		counter++;
		licztot++;
		scale=(smax-smin)/2;
		dgauss40(func,param,smin,smax,accum20,accum10);
		val20=Kahn_total(accum20);
		val10=Kahn_total(accum10);
		if(!std::isfinite(val10*val20)){
			goto Error;
		}
		
		error_ratio=fabs( (val20-val10)/(eps*val20) );
		if((error_ratio<1 ) || fabs(val20-val10)< epsabs){
			
			//total=Kahn_Sum(total,val20,accum_total,3);
			accum+=accum20;
			licz++;
			counter=0;
			
			if(smax==max){
				Kahn_free(accum);
				Kahn_free(accum10);
				Kahn_free(accum20);
				return(Kahn_total(accum));
			}else{
				smin=smax;
				increase=(4*scale);
				smax=((max-(smin+increase)<(increase/2))?(max):(smin+increase));
				
			}
		}else{
			//smax=smin+(scale/4);
			smax=smin+(scale/(1+pow(2*counter,2)*pow(error_ratio, 0.25 )));
		}
		if(((max-min)-(smax-smin))==(max-min)){
			printf("DGAUSS::division exceeds limitation. in the domain [%.3e, %.3e] of [%.3e, %.3e] after %d / %d \n",smin,smax,min,max, licz,licztot);						
			printf("valfull= %.3e , valhalf= %.3e  diff=%.3e\n",val20,val10,val20-val10);
			//getchar();
			Kahn_free(accum);
			Kahn_free(accum10);
			Kahn_free(accum20);
			goto Error;
		}
		

	}
	Error:	
		double x;
		const int N=20;
		const double *x1=x40,*x2=x20;
		for(int i=0;i<N;i++){
			x=(smax+smin)/2-scale*x1[N-i-1];
			//printf("f(%.3e) = %.3e\n",x,func(x,param));
			printf("%.3e\t%.3e\n",x,func(x,param));
		}
		x=(smax+smin)/2;
		printf("%.3e\t%.3e\n",x,func(x,param));
		for(int i=0;i<N;i++){
			x=(smax+smin)/2+scale*x1[i];
			printf("%.3e\t%.3e\n",x,func(x,param));
		}
		printf("\n");
		for(int i=0;i<N/2;i++){
			x=(smax+smin)/2-scale*x2[N/2-i-1];
			printf("%.3e\t%.3e\n",x,func(x,param));
		}
		x=(smax+smin)/2;
		printf("%.3e\t%.3e\n",x,func(x,param));
		for(int i=0;i<N/2;i++){
			x=(smax+smin)/2+scale*x2[i];
			printf("%.3e\t%.3e\n",x,func(x,param));
		}
		
		printf("\n");
		getchar();
		return 0;
}



//inline
 //double dgauss20(double (*func)(double*,void*), void* param, double min, double max,double  eps){
//	return(dgaussN(func,min,max,eps,20));	
//}
//inline
//double dgauss40(double (*func)(double*,void*),void* param,double min, double max,double  eps){
//	return(dgaussN(func,min,max,eps,40));
//}
//////////////////////////////////////////////////////////
//
//
/*
double func(double *X,void*){
	double x=*X;		
	return(x*x*exp(-x*x));
}

int main(){
	double* dummy;
	double val=dgauss(&func ,(void*)dummy, 0,100,1.0e-10);
	printf(" %.5e\n",val);
	return(0);
}*/

#endif
