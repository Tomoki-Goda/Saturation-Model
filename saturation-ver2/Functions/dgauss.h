#ifndef DGAUSS_H
#define DGAUSS_H
#include<math.h>
#include<stdio.h>
#include"./Kahn.h"

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




double dgauss(double (*func)(double*,void*),void* param, double min, double max,double  eps){
	double smin,smax; //section
	double x; 
	double total=0, val10=0, val20=0;
	double scale, mid;
	double prec=1.0e-14;
	double accum[3]={0}, accum_total[3]={0};
	double error=0;
	double arg1,arg2,increase;

	smax=max;
	smin=min;
	int licz=0,licztot=0;
	//double f[N/2+1];
	double f[11];
	if(fabs(min-max)<1.0e-15){
		return(0);
	}		
	while(1){
		licztot++;
		scale=(smax-smin)/2;
		if(scale<2.0e-15){
			printf("DGAUSS :: division exceeds limitation. in the domain [%.3e, %.3e] of [%.3e, %.3e] scale = %.5e\n",smin,smax,min,max,scale);
			getchar();
		}

		mid=(smax+smin)/2;

//		if(N==20){
			val10=0;
			//accum=0;
			Kahn_init(accum,3);
			for(int i=0;i<5;i++){
				x=scale*x10[i];
				arg1=mid+x;
				arg2=mid-x;

				//val10=Kahn(val10,w10[i]*( (*func)(&arg1,param)+(*func)(&arg2,param) ),&accum);
				val10=Kahn(val10,w10[i]*( (*func)(&arg1,param)+(*func)(&arg2,param) ),accum,3);

				//val10+=w10[i]*( (*func)(mid+x)+(*func)(mid-x) ) ;
			}
			//val10+=accum;
			val10=Kahn_total(val10,accum,3);
			//accum=0;
			//val10*=scale;
			val20=0;
			Kahn_init(accum,3);
			for(int i=0;i<10;i++){
				x=scale*x20[i];
				arg1=mid+x;
				arg2=mid-x;
				//val20=Kahn(val20,w20[i]*( (*func)(&arg1,param)+(*func)(&arg2,param) ),&accum);
				val20=Kahn(val20,w20[i]*( (*func)(&arg1,param)+(*func)(&arg2,param) ),accum,3);
				//val20+=w20[i]*( (*func)(mid+x)+(*func)(mid-x) ) ;
			}
			//val20+=accum;
			val20=Kahn_total(val20,accum,3);
			//accum=0;
/*		}else if(N==40){
			val10=0;
			accum=0;
			for(int i=0;i<10;i++){
				x=scale*x20[i];
				arg1=mid+x;
				arg2=mid-x;
				val10=Kahn(val10,w20[i]*( (*func)(&arg1,param)+(*func)(&arg2,param) ),&accum);
				//val10+=w20[i]*( (*func)(mid+x)+(*func)(mid-x) ) ;
			}
			val10+=accum;
			accum=0;
			//val10*=scale;
			val20=0;
			for(int i=0;i<20;i++){
				x=scale*x40[i];
				arg1=mid+x;
				arg2=mid-x;
				val20=Kahn(val20,w40[i]*( (*func)(&arg1,param)+(*func)(&arg2,param) ),&accum);
				//val20+=w40[i]*( (*func)(mid+x)+(*func)(mid-x) ) ;
			}
			val20+=accum;
		}
*/
		//val20*=scale;
		//printf("[%.5e ,  %.5e] \t %.5e \t res=%.5e  +-  %.5e, accumlator= %.3e\n",smin,smax, scale ,scale*val20,scale*fabs(val10-val20),accum_total);
		val10*=scale;
		val20*=scale;
		if(fabs(val20-val10)< eps*(1+fabs(val20) )){
			total=Kahn(total,val20,accum_total,3);
			licz++;
			
			if(fabs(smax-max)<1.0e-15){
				//printf("Efficiency %d/%d = %.3f\n",16*licz,16*licztot,((double)licz)/licztot);
				return(Kahn_total(total,accum_total,3 ));
			//bereak;
			}else{
				smin=smax;
				//smax=max;
				increase=(4*scale);
				smax=((max-(smin+increase)<(increase/2))?(max):(smin+increase));
			}
		}else{
			//smax=mid;
			smax=smin+(scale/2);
		}
	}
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
