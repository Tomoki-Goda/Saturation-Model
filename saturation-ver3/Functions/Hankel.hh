
class Hankel_aF{
		Laplacian_Sigma sigma;
		//Interpolation
		double max_prev=0;
		//Dipole_Gluon aF;
		int kt2_npts,x_npts;
		gsl_interp_accel *x_accel_ptr, *kt2_accel_ptr;
		gsl_spline2d *  spline_ptr;
		double *r_array,*sigma_array,*kt2_array,*x_array,*aF_array;
		//double kt2min=1.0e-15,kt2max=-1;
		//Hankel
		const double rmax=R_MAX;
		//const int n=1000;
		gsl_dht* trans;
		int approximate(){
			std::vector<double> par(2,0);
#if IBP==0
			gsl_dht_init(trans,0,rmax);
#else
			gsl_dht_init(trans,1,rmax);
#endif
			for(int i=0;i<kt2_npts;i++){
				r_array[i]=gsl_dht_x_sample(trans,i);
				kt2_array[i]=gsl_dht_k_sample(trans, i);
				printf("%.3e\n",r_array[i]);
			}
			//printf("start\n");
			for(int j=0;j<x_npts;j++){
				x_array[j]=pow(10,-10+10*((double)j/(x_npts+1) ));
				sigma.set_kinem(x_array[j]);
				//printf("%d\n",j);
				for(int i=0;i<kt2_npts;i++){
					par[1]=kt2_array[i];
					//printf("sigma(x= %.3e, r= %.3e)\n",x_array[j],r_array[i]);
					sigma_array[i]=sigma(r_array[i],par);
					//printf("sigma(x= %.3e, r= %.3e)= %.3e\n",x_array[j],r_array[i],sigma_array[i]);
				}
				//for(int i=0;i<kt2_npts;i++){
				//	printf("%.3e\t%.3e\n",kt2_array[i],aF_array[j*kt2_npts+i] );
				//}
				gsl_dht_apply(trans,sigma_array,aF_array+j*kt2_npts );
			}
			
			gsl_spline2d_init (spline_ptr,kt2_array, x_array, aF_array, kt2_npts, x_npts);
			//printf("Hankel done\n");
			//getchar();
			return 0;
		}
	public:
		Hankel_aF(){
		}
		~Hankel_aF(){
			free_approx();
		}
		int export_grid(FILE* file){
			double val;
			for(int i=0;i<x_npts;i++){
			for(int j=0;j<kt2_npts;j++){
				val=gsl_spline2d_eval_extrap(spline_ptr, sqrt(kt2_array[j]), x_array[i], kt2_accel_ptr, x_accel_ptr);
				fprintf(file,"%.5e\t%.5e\t%.5e\n",x_array[i],kt2_array[j],val);
			}
			}	
			return 0;	
		}
		void free_approx(){
			gsl_spline2d_free (spline_ptr);
			gsl_interp_accel_free (x_accel_ptr);
			gsl_interp_accel_free (kt2_accel_ptr);
			free(kt2_array);
			free(x_array);
			free(sigma_array);
			free(r_array);
			free(aF_array);
			gsl_dht_free(trans);
		}
		void set_max(double kt2max){
			//this->kt2max=kt2max;
			approximate();
		}
		void init(const int npts1,const int npts2,const int npts3,const double *par ){
			
			x_npts=npts1;
			kt2_npts=npts2*100;
			
			trans=gsl_dht_alloc(kt2_npts);
			r_array=(double*)malloc(kt2_npts*sizeof(double));
			sigma_array=(double*)malloc(kt2_npts*sizeof(double));
			x_array=(double*)malloc(x_npts*sizeof(double));
			kt2_array=(double*)malloc(kt2_npts*sizeof(double));
			aF_array=(double*)malloc(x_npts*kt2_npts*sizeof(double));
			x_accel_ptr = gsl_interp_accel_alloc ();
			kt2_accel_ptr = gsl_interp_accel_alloc ();
			spline_ptr = gsl_spline2d_alloc(gsl_interp2d_bicubic,kt2_npts, x_npts); 
		//	aF.init(N_APPROX/2+50,par);
			sigma.init(npts3,par,'l');
		}
		double operator()(const double x,const double kt2,const double mu2)const{			
			double val = 0;
			val=gsl_spline2d_eval_extrap(spline_ptr, sqrt(kt2), x, kt2_accel_ptr, x_accel_ptr);
		//	printf("aF(%.3e, %.3e)= %.3e\n",x,kt2,val);
		//	getchar();
			return(3.0/(8*PI*PI)*val);
		}
		
};
