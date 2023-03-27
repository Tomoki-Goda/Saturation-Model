//class Approx_aF;
//typedef struct {int i, j; Approx_aF* ptr; } parallel_arg;

/////////////////////////////////////////////////////
// Approx_aF<gluon> af //gluon is Dipole_Gluon or Gluon_GBW
// 
template<typename GLU >class Approx_aF{
	private:
		GLU *aF;
		double max_prev=0;
		
		int kt2_npts,x_npts;
		gsl_interp_accel *x_accel_ptr, *kt2_accel_ptr;
		gsl_spline2d *  spline_ptr;
		double *kt2_array=NULL,*x_array=NULL,*aF_array=NULL;
		double mu02=0;
		double sigma_0=0;
		double kt2min=1.0e-9,kt2max=-1;
		
		int alloc_flag=0;
		
		void free_approx(){
			if(alloc_flag!=0){
				gsl_spline2d_free (spline_ptr);
				gsl_interp_accel_free (x_accel_ptr);
				gsl_interp_accel_free (kt2_accel_ptr);
				free(kt2_array);
				free(x_array);
				free(aF_array);
				alloc_flag=0;
			}else{
				printf("Approx_aF cannot free\n");
			}
			
		}
		void alloc(int x_npts,int kt2_npts){
			if(alloc_flag!=1){
				x_array=(double*)malloc(x_npts*sizeof(double));
				kt2_array=(double*)malloc(kt2_npts*sizeof(double));
				aF_array=(double*)malloc(x_npts*kt2_npts*sizeof(double));
				x_accel_ptr = gsl_interp_accel_alloc ();
				kt2_accel_ptr = gsl_interp_accel_alloc ();
				//spline_ptr = gsl_spline2d_alloc(gsl_interp2d_bicubic,kt2_npts, x_npts);
				spline_ptr = gsl_spline2d_alloc(gsl_interp2d_bilinear,kt2_npts, x_npts);
				alloc_flag=1;
			}else{
				printf("Approx_aF cannot allocate\n");
			}
		}
		int approximate(const double kt2max){
			this->kt2max=kt2max;
			clock_t time=clock();
			std::chrono::system_clock walltime;
			std::chrono::time_point start= walltime.now();
			//double kt2,x;

			for (int j = 0; j < x_npts; ++j){
				double x=pow(10,-8+8*((double)j)/(x_npts-1));
				x_array[j] = x;
//#if SIGMA_APPROX==-2||SIGMA_APPROX==1
				aF->set_x(x);	
///#endif
				
#pragma omp parallel 
{
#pragma omp for schedule(dynamic)
				for(int i=0;i<kt2_npts;++i){
					//printf("%d\t/%d",i+1,kt2_npts);
					//fflush(stdout);
					double kt2=((double)i)/(kt2_npts-1);
					kt2=kt2min*pow(4*kt2max/kt2min,kt2)/2;
					kt2_array[i] = kt2;
					aF_array[i+ j*kt2_npts] = (*aF)(x,kt2,0);
					//printf("\033[2K\r");
					//fflush(stdout);
				}
}
				//printf("\033[2K\r");
				if(j!=0){
					printf("\033[1A\033[2K\r");
				}
				printf("approxed x=%.2e\n", x);
				fflush(stdout);
			}
			//printf("\033[1A\033[2K Grid done\n");
			printf("\033[1A\033[2K\r");
			gsl_spline2d_init (spline_ptr,kt2_array, x_array, aF_array, kt2_npts, x_npts);
			//}
			time-=clock();
			std::chrono::duration<double> wtime=walltime.now()-start;
			std::cout<< -((double)time/CLOCKS_PER_SEC)<< " CPU seconds " <<wtime.count()<<" seconds to approx"<<std::endl;
		//	printf("%.2e sec to approx\n",-((double)time/CLOCKS_PER_SEC) );
			return(0);
		}
/*
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
			printf("\033[1A\033[2K\r");
			gsl_spline2d_init (spline_ptr,kt2_array, x_array, aF_array, kt2_npts, x_npts);
			time-=clock();
			std::chrono::duration<double> wtime=walltime.now()-start;
			std::cout<< -((double)time/CLOCKS_PER_SEC)<< " CPU seconds " <<wtime.count()<<" seconds to approx"<<std::endl;
			//printf("%.2e sec to approx\n",-((double)time/CLOCKS_PER_SEC) );
			return(0);
		}
			
*/		
	public:
		double saturation(double x,double kt2_start){
			double val;
			double kt2=kt2_start;
			double diff=1.0e-1;
			double valprev=0;
			int flag=0;
			int counter=0;
			double grad=0;
			while(counter++<200){
				val=gsl_spline2d_eval_deriv_x(spline_ptr,kt2, x,kt2_accel_ptr, x_accel_ptr);
				if(fabs(val)<1.0e-5){
					return kt2;
				}else if(val>0&&flag==0){
					kt2+=diff;
					diff*=1.5;
				}else{
					flag=1;
					grad=(val-valprev)/diff;
					diff=-val/grad;
					kt2+=diff;
				}
				valprev=val;
			}
			printf("FAILED to find peak derivative is %.3e at %.3e\n",val,kt2);
			return 0;
		}
		int export_grid(FILE*file)const{
			for(int j=0;j< x_npts;j++){
				for(int i=0;i<kt2_npts;i++){
					//fprintf(file ,"%.10e\t%.10e\t%.10e\n",x_array[j],kt2_array[i],aF_array[i+j*kt2_npts]/sigma_0);
					fprintf(file ,"%.10e\t%.10e\t%.10e\n",log10(x_array[j]),log10(kt2_array[i]), aF_array[i+j*kt2_npts] );
				}
			}	
			return 0;
		}
	/*	Approx_aF(const Approx_aF& rhs){
			aF=rhs.aF;
			x_npts=rhs.x_npts;
			kt2_npts=rhs.kt2_npts;
			alloc(x_npts,kt2_npts);
			sigma_0=rhs.sigma_0;
			mu02=rhs.mu02;
			
			if(rhs.alloc_flag!=0){
				for(int i=0;i<x_npts;++i){
					x_array[i]=rhs.x_array[i];
					
					for(int j=0;j<kt2_npts;++j){
							aF_array[j+kt2_npts*i]=	rhs.aF_array[j+kt2_npts*i];			
					}
				}
				for(int j=0;j<kt2_npts;++j){
					kt2_array[j]=rhs.kt2_array[j];			
				}
				gsl_spline2d_init (spline_ptr,kt2_array, x_array, aF_array, kt2_npts, x_npts);
			}
			
			
		}
	 */
		Approx_aF(GLU& g){
			aF=&g;
		}
		~Approx_aF(){
			free_approx();
		}
		void set_max(double kt2max){
			this->kt2max=kt2max;
			approximate(kt2max);
			//approximate_thread(kt2max);
		}
		void init(const int npts1, const int npts2, const double * const &par){
			//aF=&glu;
			x_npts=npts1;
			kt2_npts=npts2;
			alloc(x_npts,kt2_npts);
			
			sigma_0=par[0];
#if MU02==0
			mu02 = par[3];
#else 
			mu02 = MU02;
#endif
		}
		double operator()(const double x,const double kt2,const double mu2)const{
			if(kt2>kt2max){printf("kt2 too large kt2max= %.1e kt2= %.1e diff=%.1e\n",kt2max,kt2,kt2max-kt2);}
			if(kt2<kt2min){printf("kt2 too small kt2min= %.1e kt2= %.1e diff=%.1e\n",kt2min,kt2,kt2min-kt2);}
			if(x>1){printf("x too large xmax= %.1e x= %.1e diff=%.1e\n",1.0,x,1.0-x);}
			if(x<1.0e-8){printf("x too small xmin= %.1e x= %.1e diff=%.1e\n",1.0e-8,x,x-1.0e-8);}			
			double val = 0;
			val=gsl_spline2d_eval(spline_ptr,kt2, x,kt2_accel_ptr, x_accel_ptr);
#if ALPHA_RUN==1
			val*=alpha(mu2+mu02)/0.2;
#endif
			return(val);
		}
};
