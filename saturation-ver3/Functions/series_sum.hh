
class Series_Sum{
	//typedef struct  Hankel;
	double eps(int n, int m, const double *list){ 
		//shanks epsilon series convergence algorithm??
		// I don't know the name but Seki-Aitken's generalization.
		if(m == -1){
			return(0);
		};
		if(m == 0){
			return(list[cycle(n)]);
		};
			return(eps(n + 1, m - 2,list) + 1/(eps(n + 1, m - 1,list) - eps(n, m - 1,list)));
	}
	int N;
	double *a;
	double *s;
	inline int cycle(int i){
		return(i-(N+2)*(i/(N+2)));
	}
	int posit=0;
	int conv=0;
	int x0,x1, x2;
	Kahn sum;
	public:
		Series_Sum(int n){
			N=2*n;
			posit=0;
			conv=0;
			a=(double*)calloc(N+2,sizeof(double));
			s=(double*)calloc(N+2,sizeof(double));
			sum=Kahn_init(3);
		}
		~Series_Sum(){
			Kahn_free(sum);
			free(a);
			free(s);
		}
		
		
		int append(double val){
			double ratio1=0,ratio2=0;
			x0=cycle(posit-2);
			x1=cycle(posit-1);
			x2=cycle(posit);
			
			
			a[x2]=val;
			sum+=val;
			s[x2]=Kahn_total(sum);
			
			if(posit>2){
				//try d'Alembert test to see it is converging.
				ratio1=a[x2]/a[x1];
				ratio2=fabs((a[x2]*a[x0])/pow(a[x1],2));
				//if ratio2 is close to 1 d'Alembert condition can be considered since now ratio is close to its limit.
				
				if(fabs(1-ratio2)<1.0e-2&& fabs(ratio1)<1 &&ratio1>0){
					conv++;//point! keep going. try a few more times.
				}else{
					conv=0;//reset
				}
				//printf("val=%.3e prev=%.3e prev2=%.3e 1-ratio1= %.2e 1-ratio2= %.2e conv=%d\n",a[x2],a[x1],a[x0],1-ratio1, 1-ratio2,conv );
				//getchar();
			}
			
			++posit;
			return(conv);
		}
		
		double extrapolate(int conv){
			double extrap=0;
			if(conv>N-1){
				extrap=eps(posit-1-N,N,s);//Now, if it is a convergeing series, we'll use Shanks accel. 
				printf("%.3e+- %.3e vs %.3e +- %.3e\n",extrap ,fabs(extrap-eps(posit-1-N-1,N,s)), s[x2],fabs(a[x2])); 
			}else{
				printf("not sufficient convergence\n");
			}
			return extrap;
		}

};















