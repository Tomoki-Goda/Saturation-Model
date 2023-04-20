#include"control.h"
#include"interpolation-gluon.hh"


void Approx_aF::free_approx(){
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
void Approx_aF::alloc(int x_npts,int kt2_npts){
	if(alloc_flag!=1){
		x_array=(double*)malloc(x_npts*sizeof(double));
		kt2_array=(double*)malloc(kt2_npts*sizeof(double));
		aF_array=(double*)malloc(x_npts*kt2_npts*sizeof(double));
		x_accel_ptr = gsl_interp_accel_alloc ();
		kt2_accel_ptr = gsl_interp_accel_alloc ();
		spline_ptr = gsl_spline2d_alloc(gsl_interp2d_bicubic,kt2_npts, x_npts);
		//spline_ptr = gsl_spline2d_alloc(gsl_interp2d_bilinear,kt2_npts, x_npts);
		alloc_flag=1;
	}else{
		printf("Approx_aF cannot allocate\n");
	}
}

//#if SUDAKOV>=1
//		int approximate(const double kt2max,const double mu2){
//#else
int Approx_aF::approximate(const double kt2max){

#if SUDAKOV>=1
	const double &mu2=*(this->mu2);
#else
	const double mu2=0;
#endif
	this->kt2max=kt2max;
	clock_t time=clock();
	std::chrono::system_clock walltime;
	std::chrono::time_point start= walltime.now();
	//double kt2,x;

	for (int j = 0; j < x_npts; ++j){
		//double x=pow(10,-8+8*((double)j)/(x_npts-1));
		double x=xmin*pow(X_MAX/xmin,((double)j)/(x_npts-1));
		
		x_array[j] = x;
#if GLUON_APPROX==1
		aF->set_x(x);	
#endif
		
#pragma omp parallel 
{
#pragma omp for schedule(dynamic)
		for(int i=0;i<kt2_npts;++i){

			double kt2=((double)i)/(kt2_npts-1);
			kt2=kt2min*pow(4*kt2max/kt2min,kt2)/2;
			kt2_array[i] = kt2;

			aF_array[i+ j*kt2_npts] = (*aF)(x,kt2,mu2);

		}
}
		//printf("\033[2K\r");
		if(j!=0){
			printf("\033[1A\033[2K\r");
		}
		printf("approxed x=%.2e mu2=%.2e\n", x,mu2);
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

double Approx_aF::saturation(double x,double kt2_start){
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

int  Approx_aF::export_grid(FILE*file)const{
	for(int j=0;j< x_npts;j++){
		for(int i=0;i<kt2_npts;i++){
			//fprintf(file ,"%.10e\t%.10e\t%.10e\n",x_array[j],kt2_array[i],aF_array[i+j*kt2_npts]/sigma_0);
#if SUDAKOV>=1
			fprintf(file ,"%.10e\t%.10e\t%.10e\t%.10e\n",log(x_array[j]),log(kt2_array[i]),log(*mu2), aF_array[i+j*kt2_npts] );
#else
			fprintf(file ,"%.10e\t%.10e\t%.10e\n",log(x_array[j]),log(kt2_array[i]), aF_array[i+j*kt2_npts] );
#endif					
		}
	}	
	return 0;
}


#if SUDAKOV>=1
void Approx_aF::set_max(double kt2max,const double& mu2){
	this->mu2=&mu2;
#else
void Approx_aF::set_max(double kt2max){
#endif
	this->kt2max=kt2max;
	approximate(kt2max);
	//approximate_thread(kt2max);
}

void Approx_aF::init(const int npts1, const int npts2, const double * const &par){
	//aF=&glu;
	x_npts=npts1;
	kt2_npts=npts2;
	alloc(x_npts,kt2_npts);
	
	sigma_0=par[0];
#if ALPHA_RUN==1
#if MU02==0
	mu02 = par[3];
#else 
	mu02 = MU02;
#endif
#endif
}

double Approx_aF::operator()(const double x,const double kt2,const double mu2)const{
#if SUDAKOV>=1
	if(*(this->mu2)!=mu2){
		printf("Approx_aF:: mu2 doesn't match. internal = %.3e, mu2=%.3e\n",*(this->mu2),mu2);
	}
#endif
	if(kt2>kt2max){printf("Approx_aF:: kt2 too large kt2max= %.1e kt2= %.1e diff=%.1e\n",kt2max,kt2,kt2max-kt2);}
	if(kt2<kt2min){printf("Approx_aF:: kt2 too small kt2min= %.1e kt2= %.1e diff=%.1e\n",kt2min,kt2,kt2min-kt2);}
	if(x>X_MAX){printf("Approx_aF:: x too large xmax= %.1e x= %.1e diff=%.1e\n",X_MAX,x,X_MAX-x);}
	if(x<X_MIN){printf("Approx_aF:: x too small xmin= %.1e x= %.1e diff=%.1e\n",1.0e-8,x,x-1.0e-8);}			
	double val = 0;
	val=gsl_spline2d_eval(spline_ptr,kt2, x,kt2_accel_ptr, x_accel_ptr);
#if ALPHA_RUN==1
	val*=alpha(mu2+mu02)/0.2;
#endif
	return(val);
}
