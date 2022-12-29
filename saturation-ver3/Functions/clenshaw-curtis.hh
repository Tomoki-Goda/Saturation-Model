//c++ version of clenshaw-curtis integration function
#ifndef CLENSHAW_LOADED
#define CLENSHAW_LOADED 1

#include<cmath>
#include<iostream>
#ifndef PI
#define PI 3.141592653589793238462643383279502884197
#endif
#include<quadmath.h>
#include"./kahnsum.hh"
#include"./CC_FIX.hh"
#include"./dgauss.h"
//#include"./Kahn.h"
//#include"./clenshaw.h"



class Clenshaw_Curtis:public Clenshaw_Curtis_FIX{
	using Clenshaw_Curtis_FIX::Clenshaw_Curtis_FIX;
	private:
		//
		double  integrate_cc(double(*func)(double*,void*),void* args,double min,double max,int sect, double eps, double epsabs,int *error)const{
			static int warn_flag=0;
			double smin,smax,dummy;
			double **par=(double**)args;
			double scale;//,mid;
			double valfull,valhalf;
			
			Kahn total(2);
			total.zero();
			if(min>max){
				printf("max < min \n");
				smin=max;
				max=min;
				min=smin;
			}

			smin=min;
			smax=min + (max-min)/sect;

			unsigned int licz=0;
			unsigned int totlicz=0;
			double efficiency,increase;
			if(min==max){
				return 0;
			}
			
			int counter=0;
			while(true){
				counter++;
				totlicz++;
				if(totlicz==50){
					printf("Ran %d times.  Reconsider strategy \n",totlicz);
				}
				scale=(smax-smin)/2;
				
			
				integral16(func,args,smin,smax,&valfull,&valhalf);

				if(Print==1){
					printf("%.5e , %.5e in the domain [%.3e, %.3e] of [%.3e, %.3e] \n",valfull,valhalf,smin,smax,min,max);
				}
				//if(fabs(valfull-valhalf)<eps*(1.0e-3+fabs(valfull)) ){//originally 1 but now 1.0e-3 no reason but to go further in the absolute accuracy
				//if(( fabs(valfull-valhalf)<eps*(fabs(valfull)) ) || (  fabs(2*valfull*(max-min)/sqrt(scale))<epsabs )){//Need improvement
				if((max-min)<1.0e-10&&valfull<eps){
					return(valfull);
				}
	
				if(( fabs(valfull-valhalf)<eps*(fabs(valfull)) ) || (  fabs(valfull-valhalf)<epsabs )||(counter>=MAX_RECURSION)){//Need improvement
					if((!(( fabs(valfull-valhalf)<eps*(fabs(valfull)) ) || (  fabs(valfull-valhalf)<epsabs ))) && counter>=MAX_RECURSION&&warn_flag<3){
							warn_flag++;
							std::cout<<name<<"\t";
							printf("MAX_RECURSION:: evaluated %d times. Increase MAX_RECURSION\n",MAX_RECURSION );
							printf("[%.3e, %.3e] of [%.3e, %.3e] after %d / %d \n",smin,smax,min,max, licz,totlicz);
							printf("valfull= %.10e , valhalf= %.10e  diff=%.3e\n",valfull,valhalf,valfull-valhalf);
							//dummy=(smax+smin)/2;
							//printf("\nf(%.3e )=%.3e f(%.3e )=%.3e f(%.3e )=%.3e\n", smin,(*func)(&smin,args),dummy,(*func)(&dummy,args),smax ,(*func)(&smax,args));
							for(int i=0;i<(N+1);i++){
								dummy=-cos(i*PI/N)*(smax-smin)+(smax+smin);
								dummy/=2;
								//printf("%d: f(%.3e )=%.3e\n",i,dummy,(*func)(&dummy,args));
								printf("%.5e\t%.5e\n",dummy,(*func)(&dummy,args));
							}
							
							//printf("beta=%.3e kappa=%.3e kt=%.3e %.3e %.3e %.3e\n\n", par[0][0],par[0][1],par[0][2],par[0][3],par[0][4],par[0][5] );
							*error=1;
							return(0 );
							if(warn_flag==2){
								printf("\033[0;31m Warning message suppressed.\033[m\n");
							}
							printf("\n");					
					}
					total+=valfull;
					licz++;
					
					if(smax==max){
						efficiency=((double)licz)/totlicz;
						if(Print==1||efficiency<0.5){
							//DIV*=2;
							std::cout<<name<<" ";
							printf("Clenshaw_Curtis:: Efficiency %d/%d=%f\n",N*licz,N*totlicz,((double)licz)/totlicz);
						}
						return(total.total());
					}
					increase=2*(DIV*counter)*scale;//twice the current section size to start.
					smin=smax;
					smax=( ( (max-(increase+smin))<increase/2)?(max):(increase+smin) );//but only if remaining section is not too small.(it is wasteful to compute small remnant of section at the end...
					counter=0;
				}else{  
					smax=smin+(scale/(DIV*counter) );
				}
				if((smax-smin)<1.0e-15){
					std::cout<<name<<" ";
					printf("Clenshaw_Curtis::division exceeds limitation. in the domain [%.3e, %.3e] of [%.3e, %.3e] after %d / %d \n",smin,smax,min,max, licz,totlicz);
					printf("valfull= %.10e , valhalf= %.10e \n",valfull,valhalf);
					getchar();
				}
			}
		}
	
	//double  integrate_cc(double(*func)(double*,void*),void* args,double min,double max,int sect, double eps)const{
	//	return (integrate_cc(func, args, min, max, sect,  eps,1,));
	//}
public:
	int ERROR=0;
	int MAX_RECURSION=5;
	std::string name="";
	
 	double	integrate(double(*func)(double*,void*),void* args,double min,double max, int sect,double eps,double epsabs){
		return( integrate_cc(func, args, min, max,sect,  eps, epsabs,&ERROR));
	}
	double  integrate(double(*func)(double*,void*),void* args,double min,double max,int sect, double eps){
		return (integrate_cc(func, args, min, max, sect,  eps, eps,&ERROR));
	}
	double	integrate(double(*func)(double*,void*),void* args,double min,double max, double eps,double epsabs){
		return( integrate_cc(func, args, min, max,1,  eps, epsabs,&ERROR));
	}
	double  integrate(double(*func)(double*,void*),void* args,double min,double max, double eps){
		return (integrate_cc(func, args, min, max, 1,  eps, eps,&ERROR));
	}

	
};
/*
 double func(double *X){
	double x=*X;
	double val=(1-50*x*x)*exp(-x*x);
	return(val);
}
int main(){
	Clenshaw_Curtis integration(16);
	double res=integration.integrate(&func,0,100,1.0e-10);
	std::cout<<res<<std::endl;
 	return(0);
	
}*/
#endif
