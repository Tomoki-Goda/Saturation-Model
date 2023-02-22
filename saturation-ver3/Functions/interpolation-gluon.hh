class Approx_aF;
typedef struct {int i, j; Approx_aF* ptr; } parallel_arg;

class Approx_aF{
	private:
		Dipole_Gluon aF;
		double max_prev=0;
		
		int kt2_npts,x_npts;
		gsl_interp_accel *x_accel_ptr, *kt2_accel_ptr;
		gsl_spline2d *  spline_ptr;
		double *kt2_array,*x_array,*aF_array;
		double mu02=0;
		
		double kt2min=1.0e-15,kt2max=-1;
		
		void free_approx(){
			gsl_spline2d_free (spline_ptr);
			gsl_interp_accel_free (x_accel_ptr);
			gsl_interp_accel_free (kt2_accel_ptr);
			free(kt2_array);
			free(x_array);
			free(aF_array);
		}
		int approximate(const double kt2max){
			this->kt2max=kt2max;
			clock_t time=clock();
			std::chrono::system_clock walltime;
			std::chrono::time_point start= walltime.now();
			double kt2,x;
			for (int j = 0; j < x_npts; ++j){
				x=pow(10,-8+8*((double)j)/(x_npts-1));
				x_array[j] = x;
				aF.set_x(x);	
				if(j!=0){
					printf("\033[1A\033[2K\r");
				}
				printf("[ ");
				for(int i=0;i<kt2_npts;++i){
					kt2=((double)i)/(kt2_npts-1);
					kt2=kt2min*pow(4*kt2max/kt2min,kt2)/2;
					kt2_array[i] = kt2;
					aF_array[i+ j*kt2_npts] = aF(kt2,0);
					if((i/4)*4==i){
						printf("=");
					}
				}
				printf("\033[2K\r");
				printf(" approxed x=%.2e\n", x);
			}
			//printf("\033[1A\033[2K Grid done\n");
				
			gsl_spline2d_init (spline_ptr,kt2_array, x_array, aF_array, kt2_npts, x_npts);
			//}
			time-=clock();
			std::chrono::duration<double> wtime=walltime.now()-start;
			std::cout<< -((double)time/CLOCKS_PER_SEC)<< " CPU seconds " <<wtime.count()<<" seconds to approx"<<std::endl;
		//	printf("%.2e sec to approx\n",-((double)time/CLOCKS_PER_SEC) );
			return(0);
		}
		static void* compute(void* par){
			//printf("func\n");
			parallel_arg* param=(parallel_arg*)par;
			int i=param->i,j=param->j;
			Approx_aF* to=param->ptr;
			double kt2=((double)i)/(to->kt2_npts-1);
			kt2=to->kt2min*pow(4*(to->kt2max)/(to->kt2min),kt2)/2;
			(to->kt2_array)[i] = kt2;
			(to->aF_array)[i+ j*(to->kt2_npts)] = (to->aF)(kt2,0);
			return NULL;
		}
		int approximate_thread(const double kt2max){
			this->kt2max=kt2max;
			clock_t time=clock();
			std::chrono::system_clock walltime;
			std::chrono::time_point start= walltime.now();
			double kt2,x;
			pthread_t thread[kt2_npts];
			parallel_arg args[kt2_npts];
			int i1;//,i2,i3,i4;
			for (int j = 0; j < x_npts; ++j){
				x=pow(10,-8+8*((double)j)/(x_npts-1));
				x_array[j] = x;
				aF.set_x(x);	
				if(j!=0){
					printf("\033[1A\033[2K\r");
				}
				printf("[ ");
				for(int i=0;i<kt2_npts;++i){
					args[i].i=i;
					args[i].j=j;
					args[i].ptr=this;
					i1=pthread_create(thread+i,NULL,compute,(void*)(args+i) );
				}					
				for(int i=0;i<kt2_npts;++i){
					pthread_join(thread[i],NULL);
				}
				printf("\033[2K\r");
				printf(" approxed x=%.2e\n", x);
			}
			//printf("\033[1A\033[2K Grid done\n");
			gsl_spline2d_init (spline_ptr,kt2_array, x_array, aF_array, kt2_npts, x_npts);
			time-=clock();
			std::chrono::duration<double> wtime=walltime.now()-start;
			std::cout<< -((double)time/CLOCKS_PER_SEC)<< " CPU seconds " <<wtime.count()<<" seconds to approx"<<std::endl;
			printf("%.2e sec to approx\n",-((double)time/CLOCKS_PER_SEC) );
			return(0);
		}

	public:
		int export_grid(FILE*file)const{
			for(int j=0;j< x_npts;j++){
				for(int i=0;i<kt2_npts;i++){
					fprintf(file ,"%.10e\t%.10e\t%.10e\n",x_array[j],kt2_array[i],aF_array[i+j*kt2_npts]);
				}
			}	
			return 0;
		}

		Approx_aF(){
		}
		~Approx_aF(){
			free_approx();
		}
		void set_max(double kt2max){
			this->kt2max=kt2max;
			//approximate(kt2max);
			approximate_thread(kt2max);
		}
		void init(const int npts1, const int npts2, const int npts3, const double *par ){
			x_npts=npts1;
			kt2_npts=npts2;
			
			x_array=(double*)malloc(x_npts*sizeof(double));
			kt2_array=(double*)malloc(kt2_npts*sizeof(double));
			//kt2_array=(double*)mmap(NULL,kt2_npts*sizeof(double),PROT_READ|PROT_WRITE,MAP_SHARED|MAP_ANONYMOUS,-1,0);
			aF_array=(double*)malloc(x_npts*kt2_npts*sizeof(double));
			//aF_array=(double*)mmap(NULL,x_npts*kt2_npts*sizeof(double),PROT_READ|PROT_WRITE,MAP_SHARED|MAP_ANONYMOUS,-1,0);
			x_accel_ptr = gsl_interp_accel_alloc ();
			kt2_accel_ptr = gsl_interp_accel_alloc ();
			spline_ptr = gsl_spline2d_alloc(gsl_interp2d_bicubic,kt2_npts, x_npts); 
			aF.init(npts3,par);
#if MU02==0
			mu02 = par[3];
#else 
			mu02 = MU02;
#endif
		}
		double operator()(const double x,const double kt2,const double mu2)const{			
			double val = 0;
			
			val=gsl_spline2d_eval(spline_ptr,kt2, x,kt2_accel_ptr, x_accel_ptr);
#if ALPHA_RUN==1
			val*=alpha(mu2+mu02)/0.2;
			//printf("%.2e %.2e\n",mu2+mu02,alpha(mu2+mu02));
#endif
			return(val);
		}
};
