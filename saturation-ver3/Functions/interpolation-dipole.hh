

class Laplacian_Sigma;
typedef struct {int j; Laplacian_Sigma *ptr;} sigmaopt;

//FOR APPROXIMATION AND DERIVATIVES		
//class Laplacian_Sigma:public Sigma{

//pthread_mutex_t mut;
class Laplacian_Sigma{
	private:
		Sigma sigma;
		double max=R_MAX, min=R_MIN; 
		double *r_array=NULL,*sigma_array=NULL;
		gsl_interp_accel *  r_accel_ptr;
		gsl_spline *  spline_ptr;
		char mode='l';//l or s
		int r_npts=0;
		int counter=0;
		double x=0;
		double sigma_0=0;
		
		void free_approx(){
			gsl_spline_free (spline_ptr);
			gsl_interp_accel_free (r_accel_ptr);
			free(r_array);
			free(sigma_array);
		}
		void allocate(int npts1){
			r_npts=npts1;
			r_array=(double*)calloc(r_npts,sizeof(double));
			sigma_array=(double*)calloc(r_npts,sizeof(double));
			r_accel_ptr = gsl_interp_accel_alloc ();
			spline_ptr = gsl_spline_alloc(gsl_interp_cspline, r_npts); // cubic spline
		}
		static void* compute(void*opt){
			sigmaopt* param=(sigmaopt*)opt;
			Laplacian_Sigma *sigmaptr=param->ptr;
			const int j=param->j;
			(sigmaptr->sigma_array)[j] = (sigmaptr->sigma)((sigmaptr->r_array)[j]);
			if(!std::isfinite((sigmaptr->sigma_array)[j])){
				printf("can not approximate sigma=%.3e\n",sigmaptr->sigma_array[j] );
			}
			return opt;
		}
		int approximate_thread(const double x){
			sigma.set_kinem(x);
			pthread_t thread[r_npts];
			sigmaopt args[r_npts];
			for (int j = 0; j < r_npts; j++){
				args[j].j=j;
				args[j].ptr=this;
				pthread_create(thread+j,NULL,&compute,(void*)(&args[j]) );
			}
			for(int i=0;i<r_npts;++i){
					pthread_join(thread[i],NULL);
			}
			gsl_spline_init (spline_ptr, r_array, sigma_array, r_npts);
			//test();
			return(0);
		}
		int approximate(const double x){
			sigma.set_kinem(x);
			for (int j = 0; j < r_npts; j++){
				sigma_array[j] = sigma(r_array[j]);
				if(!isfinite(sigma_array[j])){
					printf("can not approximate sigma=%.3e\n",sigma_array[j] );
				}
			}
			gsl_spline_init (spline_ptr, r_array, sigma_array, r_npts);
			return(0);
		}
		int test(){
			double diff[2*r_npts-1];
			double sigma_array_2[2*r_npts-1];
			for (int j = 0; j < r_npts-1; j++){
				sigma_array_2[j] = sigma( (r_array[j]+2*r_array[j+1])/3 );
				diff[j]=sigma_array_2[j]-gsl_spline_eval(spline_ptr, (r_array[j]+2*r_array[j+1])/3 ,r_accel_ptr);
				//sigma_array_2[j] = sigma( r_array[j] );
				//diff[j]=sigma_array_2[j]-gsl_spline_eval(spline_ptr, r_array[j] ,r_accel_ptr);
				if(!isfinite(sigma_array[j])){
					printf("can not approximate sigma=%.3e\n",sigma_array[j] );
				}
				if(fabs(diff[j]/sigma_array_2[j])>1.0e-5&& abs(diff[j])>1.0e-10){
					printf("%.3e\t%.3e\t%.3e\t%.3e\n",x,r_array[j],sigma_array_2[j],diff[j]);
				}
			}printf("\n");
			//exit(1);
			//gsl_spline_init (spline_ptr, r_array, sigma_array, r_npts);
			return(0);
		}
		
	public:
		//double max=R_MAX, min=R_MIN;
		inline int set_kinem(double x){
			this->x=x;
			//approximate(x);
			approximate_thread(x);
			return 0;
		}	
		
		explicit Laplacian_Sigma(){
			r_npts=0;
			sigma_array=NULL;
			r_array=NULL;
			//printf("sigma approx\n");
		}
		~Laplacian_Sigma(){
			//printf("sigma approx end\n");
			free_approx();
			
		}
		void init(const int npts1,const double *par ,char mode){
			sigma_0=par[0];
			this->mode=mode;
			if(r_npts!=0){
				free_approx();
			}
			allocate(npts1);
			sigma.init(par);
			
			double r;
			for (int j = 0; j < r_npts; j++){
				r=((double)j)/(r_npts-1);
				r=min*pow(4.0*max/min,r)/2.0;
				r_array[j]=r;
			}
		}
		double operator()(const double rho)const{
			double var;
#if R_CHANGE_VAR==1
			const double r=rho/(1-rho);
			var=pow(1-rho,-2);
#elif R_CHANGE_VAR==0
			const double r =rho;
			var=1;
#endif		
			var=gsl_spline_eval(spline_ptr, r,r_accel_ptr);
		
			return(var);
		}
		//int export_grid(FILE* file, FILE* file2){
		int export_grid(FILE* file){
			for(int i=0;i<r_npts;i++){
				fprintf(file,"%.5e\t%.5e\t%.5e\n",x,r_array[i],sigma_array[i]/sigma_0);
			}		
		return 0;	
		}
		double operator()(const double rho, const std::vector<double> &par)const {
#if R_CHANGE_VAR==1
			const double r=rho/(1-rho);
#elif R_CHANGE_VAR==0
			const double r =rho;
#endif		
			if(r>2*R_MAX){
				printf("too large, out of range %.4e - %.4e = %.4e\n",R_MAX,r,R_MAX-r);
			}
			if(r<R_MIN/2){
				printf("too small, out of range %.4e - %.4e = %.4e , rho=%.3e,%d\n",r,R_MIN,r-R_MIN,rho,R_CHANGE_VAR);
			}
			const double kt2=par[0],x=par[1];
			double val = 0;

//#if (LAPLACIAN==1||R_FORMULA==1)
			switch(mode){
				case 'l':
#if IBP==1
					val=gsl_spline_eval_deriv(spline_ptr, r,r_accel_ptr);
					val*= sqrt(kt2)*r*std::cyl_bessel_j(1,r*sqrt(kt2));
#elif IBP==2
					val=gsl_spline_eval(spline_ptr, r,r_accel_ptr);
					val*=-kt2*r*std::cyl_bessel_j(0,r*sqrt(kt2));
#else
					val=gsl_spline_eval_deriv2(spline_ptr, r,r_accel_ptr);
					val+=gsl_spline_eval_deriv(spline_ptr, r,r_accel_ptr)/r;
					val*=r*std::cyl_bessel_j(0,r*sqrt(kt2));
#endif			
					break;
				case 's':
					val=gsl_spline_eval(spline_ptr, r,r_accel_ptr);
					val*=r*std::cyl_bessel_j(0,r*sqrt(kt2));
					break;
				default:
					printf("unknown option in laplacian sigma\n");
			}
			//printf("2: val=%.3e for r= %.3e\n",val,r);
			if(not(std::isfinite(val))){
				printf("2: val=%.3e for r= %.3e\n",val,r);
			}
			
			if(not(std::isfinite(val))){
				printf("3: val=%.3e for r= %.3e\nparameter: ",val,r);
				for(int i=0;i<par.size();i++){
					printf(" %.3e\t",par[i]);
				}printf("\n");
			}
#if (R_CHANGE_VAR==1)
			return(val/pow(1-rho,2));
#elif (HANKEL==1||R_CHANGE_VAR==0)
			return val;
#endif
		}
		
		double constant(double r , const std::vector<double> &par)const {
			//const double r=rho/(1-rho);
			const double kt=sqrt(par[0]);
			double val=0;
#if IBP==1
			val=gsl_spline_eval_deriv(spline_ptr, r,r_accel_ptr);
			val*=r*std::cyl_bessel_j(0,r*kt);
#elif IBP==2
			val=kt*std::cyl_bessel_j(1,r*kt)*(gsl_spline_eval(spline_ptr,r,r_accel_ptr));
			val+=std::cyl_bessel_j(0,r*kt)*gsl_spline_eval_deriv(spline_ptr,r,r_accel_ptr);
			val*=r;
#endif
			return(val);
		}
};





/*
double par[]={25.0, 0.3,3.0e-4};
double Q0=1;
#define R_MIN 1.0e-6
#define R_MAX 50
double INT_PREC=1.0e-5;

int main(){
	FILE* file=fopen("./test.txt","w");
	double Q2=10;
	double x=0.01;
	double val,kt2=0;
	Approx_aF gbw;
	gbw.init(100,100,par);
	gbw.set_max(1.0e+3);
	//std::vector<double> param(2,0);
	//param[0]=x;
	//param[1]=1;

	for(int i =0;i<500;i++){
		kt2=pow(10,-4+7*((double)i)/(500));
		//val=gbw(kt2,param);
		val=gbw(x,kt2);
		fprintf(file,"%.3e\t%.3e\n",kt2,val);
		//printf("%.3e\t%.3e\n",kt2,val);
	}
	
	fclose(file);
	
		
	
}
*/
