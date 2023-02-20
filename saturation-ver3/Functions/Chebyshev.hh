double func_cheb(const double rho, void*p){
	const double r=rho/(1-rho);
	Sigma* sigma=(Sigma*)p;		
	double val=(*sigma)(r);///pow(1-rho,2);
	//printf("func %.3e r=%.3e, rho=%.3e\n",val,r, rho);getchar();
	return val;
}
double func_deriv(const double rho, void*p){
	gsl_function F;
 	F.function = func_cheb;
 	F.params=p;
	//const double r=rho/(1-rho);
	double val,err;
	gsl_deriv_forward(&F,rho,rho/10, &val, &err);
	
	return val;
}

class Laplacian_Sigma_Cheb{

	private:
		
		Sigma sigma;
		gsl_cheb_series *coeffs,*deriv;
		int r_npts=0;
		double x=0;
		char mode='l';
		void free_approx(){
			gsl_cheb_free (coeffs);
			gsl_cheb_free (deriv);
		}
		
	public:
		//double max=R_MAX, min=R_MIN;
		inline int set_kinem(double x){
			approximate(x);
			return 0;
		}	
		int approximate(const double x){
			this->x=x;
			sigma.set_kinem(x);
			//printf("Approximate\n");
			//getchar();
			gsl_function F;
 			F.function = func_cheb;
 			//F.function = func_deriv;
 			F.params=(void*)&sigma;
			//gsl_cheb_init(coeffs,&F,R_MIN,R_MAX);
			gsl_cheb_init(coeffs,&F,R_MIN/(1+R_MIN),R_MAX/(1+R_MAX));
			gsl_cheb_calc_deriv(deriv,coeffs);
			//printf("End\n");
			//getchar();
			return(0);
		}
		explicit Laplacian_Sigma_Cheb(){
			r_npts=0;
			
			//printf("sigma approx\n");
		}
		~Laplacian_Sigma_Cheb(){
			//printf("sigma approx end\n");
			free_approx();
			
		}
		void init(const int npts1,const double *par ,char mode){
			this->mode=mode;
			
			if(r_npts!=0){
				gsl_cheb_free(coeffs);
				gsl_cheb_free(deriv);
			}//else{
				//printf("first\n");
			r_npts=npts1;
			coeffs=gsl_cheb_alloc(npts1);
			deriv=gsl_cheb_alloc(npts1);
			sigma.init(par);
		}
		double operator()(const double rho)const{
			//const double r=rho/(1-rho);
			return(gsl_cheb_eval(deriv,rho));///pow(1-rho,2));
		}
		double operator()(const double rho, const std::vector<double> &par)const {
			const double r=rho/(1-rho);
			
			const double kt2=par[1];
			double val = 0;

//#if (LAPLACIAN==1||R_FORMULA==1)
			switch(mode){
				case 'l':
#if (IBP>=1 && HANKEL!=1)
#if IBP==1
					val=gsl_cheb_eval(deriv,rho);
					//val=gsl_cheb_eval(coeffs,rho);
					val*= sqrt(kt2)*r*std::cyl_bessel_j(1,r*sqrt(kt2));
#elif IBP==2
					val=gsl_cheb_eval(coeffs,rho);
					//val=gsl_cheb_eval(coeffs,rho);
					val*= kt2*r*std::cyl_bessel_j(0,r*sqrt(kt2));
#endif
#elif HANKEL==1
					val=gsl_cheb_eval(deriv,r);;
#else
					val=0;
#endif			
					break;
				case 's':
					val=gsl_cheb_eval(coeffs,r);
					val*=r*std::cyl_bessel_j(0,r*sqrt(kt2));
					break;
				default:
					printf("unknown option in laplacian sigma\n");
			}
		
			if(not(std::isfinite(val))){
				printf("2: val=%.3e for r= %.3e\n",val,r);
			}
			
			if(not(std::isfinite(val))){
				printf("3: val=%.3e for r= %.3e\nparameter: ",val,r);
				for(int i=0;i<par.size();i++){
					printf(" %.3e\t",par[i]);
				}printf("\n");
			}
			return(val);
		}
		
		double constant(double rho , const std::vector<double> &par)const {
			const double r=rho/(1-rho);
			const double kt2=par[1];
			double val;
			val=gsl_cheb_eval(deriv,rho)*pow(1-rho,2);
			val*=r*std::cyl_bessel_j(0,r*sqrt(kt2));
		//	if(fabs(val)>1.0e-3){
		//		printf("x= %.2e boundar= %.2e\n",x,val);
		//	}
			return(val);
		}
};
