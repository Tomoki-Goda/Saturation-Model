
#include"control.h"
#include"gluons.hh"
#include"r-formula.hh"
#include"gluon-integrand.hh"
#include"miscellaneous.hh"
/////////////////////////////////////////////////////////////////////
//Gluon_Integrand
////////////////////////////////////////////////////////////////////
template <typename Sig>void Gluon_Integrand<Sig>::init(const double * const &par ,char mode){
	this->mode=mode;
	//sigma->init(par);
}
template <typename Sig>int Gluon_Integrand<Sig>::set_x(const double &x){
	sigma->set_x(x);
	fixx=&x;
	return 0;
}
template <typename Sig>double Gluon_Integrand<Sig>::operator()(const double rho, const std::vector<double> &par) {
#if R_CHANGE_VAR==1
	const double r=rho/(1-rho);
#elif R_CHANGE_VAR==0
	const double r =rho;
#endif		
	if(r>1.42*R_MAX){
		printf("Gluon_Integrand:: too large, out of range %.4e - %.4e = %.4e\n",R_MAX,r,R_MAX-r);
	}
	if(r<R_MIN/1.41){
		printf("Gluon_Integrand:: too small, out of range %.4e - %.4e = %.4e , rho=%.3e,%d\n",r,R_MIN,r-R_MIN,rho,R_CHANGE_VAR);
	}
	const double kt=sqrt(par[0]),x=par[1];
#if SUDAKOV>=1
	const double q2=par[2];
#endif
	double val = 0;
	if(x!=*fixx){
		printf("Gluon_Integrand:: Error: x does not match. input=%.3e internal x= %.3e diff = %.3e\n",x,*fixx, x-*fixx);
	}
	double h=r/50;
	switch(mode){
		case 'l':
#if IBP==1
			val=deriv<Sig>(*sigma,x,r,h,1);
#if NS==2
			val*=( kt*r*std::cyl_bessel_j(1,r*kt) + 2*pow(r/ns_pow,2)*std::cyl_bessel_j(0,r*kt));
#else
			val*= kt*r*std::cyl_bessel_j(1,r *kt );
#endif
#elif IBP==2
			val=(*sigma)(x,r);
			val*=-kt*kt*r*std::cyl_bessel_j(0,r*kt);
#else
			//The following line is very inefficient!!
			val=deriv<Sig>(*sigma,x,r,h,2);
			val+=deriv<Sig>(*sigma,x,r,h,1)/r;
			val*=r*std::cyl_bessel_j(0,r*kt);
#endif			
			break;
		case 's':
			val=(*sigma)(x,r);
			val*=r*std::cyl_bessel_j(0,r*kt);
			break;
		case 'w'://for weizsacker-william 
			val=(*sigma)(x,r);
			val*=std::cyl_bessel_j(0,r *kt );
#if SUDAKOV>=1		
			val/=exp(sudakov(r,q2,kt*kt));//exp(-S);	
#endif
			val/=r;
			break;
		default:
			printf("Gluon_Integrand:: unknown option in laplacian sigma\n");
	}
	//printf("2: val=%.3e for r= %.3e\n",val,r);
	if(not(std::isfinite(val))){
		printf("Gluon_Integrand:: 2: val=%.3e for r= %.3e\n",val,r);
	}
	
	if(not(std::isfinite(val))){
		printf("Gluon_Integrand:: 3: val=%.3e for r= %.3e\nparameter: ",val,r);
		for(int i=0;i<par.size();i++){
			printf(" %.3e\t",par[i]);
		}printf("\n");
	}
#if NS>=1
	val*=exp(-pow(r/ns_pow,2));
#endif
#if (R_CHANGE_VAR==1)
	return(val/pow(1-rho,2));
#elif (HANKEL==1||R_CHANGE_VAR==0)
	return val;
#endif
}

template <typename Sig>double Gluon_Integrand<Sig>::constant(double r , const std::vector<double> &par) {
#if WW==1
return 0;
#endif
	const double kt=sqrt(par[0]),x=par[1];
	double val=0;
	double h=r/50;
	if(x!=*fixx){
		printf("Gluon_Integrand:: Error: x does not match. input=%.3e internal x= %.3e diff = %.3e\n",x,*fixx, x-*fixx);
	}
#if IBP==1
	val=deriv< Sig>(*sigma,x,r,h,1);
	val*=r*std::cyl_bessel_j(0,r*kt);
#elif IBP==2
	val=kt*std::cyl_bessel_j(1,r*kt)*(*sigma)(x,r);
	val+=std::cyl_bessel_j(0,r*kt)*deriv< Sig>(*sigma,x,r,h,1);
	val*=r;
#endif
#if NS>=1
	val*=exp(-pow(r/ns_pow,2));
#endif
	return(val);
}



////////////////////////////////////////////////////////////////////
//Laplacian_Sigma
// This is kept only for the testing purpose, no-longer maintained
///////////////////////////////////////////////////////////////////
template <typename Sig>void Laplacian_Sigma<Sig>::free_approx(){
	if(alloc_flag!=0){
	gsl_spline_free (spline_ptr);
	gsl_interp_accel_free (r_accel_ptr);
	free(r_array);
	free(sigma_array);
	alloc_flag=0;
	}else{
		printf("Laplacian_Sigma cannot free\n");
	}
}
template <typename Sig>void Laplacian_Sigma<Sig>::allocate(int npts1){
	if(alloc_flag!=1){
		r_npts=npts1;
		r_array=(double*)calloc(r_npts,sizeof(double));
		sigma_array=(double*)calloc(r_npts,sizeof(double));
		r_accel_ptr = gsl_interp_accel_alloc ();
		spline_ptr = gsl_spline_alloc(gsl_interp_cspline, r_npts); // cubic spline
		//spline_ptr = gsl_spline_alloc(gsl_interp_akima, r_npts); 
		//spline_ptr = gsl_spline_alloc(gsl_interp_linear, r_npts);
		//spline_ptr = gsl_spline_alloc(gsl_interp_steffen, r_npts);
		alloc_flag=1;
	}else{
		printf("Laplacian_Sigma cannot allocate\n");
	}
}
template <typename Sig>int Laplacian_Sigma<Sig>::approximate(const double x){
#pragma omp parallel
{
#pragma omp for schedule(dynamic)
	for (int j = 0; j < r_npts; j++){
		sigma_array[j] = (*sigma)(x,r_array[j]);
		if(!isfinite(sigma_array[j])){
			printf("can not approximate sigma=%.3e\n",sigma_array[j] );
		}
	}
}
	gsl_spline_init (spline_ptr, r_array, sigma_array, r_npts);
	return(0);
}

template <typename Sig>inline int Laplacian_Sigma<Sig>::set_x(const double& x){
	fixx=&x;
	sigma->set_x(x);
	approximate(x);
	return 0;
}	

template <typename Sig>void Laplacian_Sigma<Sig>::init(const int npts1,const double * const &par ,char mode){
	sigma_0=par[0];
	this->mode=mode;
	if(r_npts!=0){
		free_approx();
	}
	allocate(npts1);
	
	double r;
	for (int j = 0; j < r_npts; j++){
		r=((double)j)/(r_npts-1);
		r=(1-cos(PI*r))/2;
		r=R_MIN*pow(2.0*R_MAX/R_MIN,r)/1.414;
		r_array[j]=r;
	}
}
template <typename Sig>double Laplacian_Sigma<Sig>::operator()(const double rho)const{
	double var;
#if R_CHANGE_VAR==1
	const double r=rho/(1-rho);
	var=pow(1-rho,-2);
#elif R_CHANGE_VAR==0
	const double r =rho;
	var=1;
#endif		
	var*=gsl_spline_eval(spline_ptr, r,r_accel_ptr);

	return(var);
}
//int export_grid(FILE* file, FILE* file2){
template <typename Sig>int Laplacian_Sigma<Sig>::export_grid(FILE* file){
	for(int i=0;i<r_npts;i++){
		fprintf(file,"%.5e\t%.5e\t%.5e\n",*fixx,r_array[i],sigma_array[i]/sigma_0);
	}		
return 0;	
}

template <typename Sig>double Laplacian_Sigma<Sig>::operator()(const double rho, const std::vector<double> &par){
#if R_CHANGE_VAR==1
	const double r=rho/(1-rho);
#elif R_CHANGE_VAR==0
	const double r =rho;
#endif		
	if(r>2*rmax){
		printf("too large, out of range %.4e - %.4e = %.4e\n",rmax,r,rmax-r);
	}
	if(r<R_MIN/2){
		printf("too small, out of range %.4e - %.4e = %.4e , rho=%.3e,%d\n",r,R_MIN,r-R_MIN,rho,R_CHANGE_VAR);
	}
	const double kt=sqrt(par[0]),x=par[1];
	if(x!=*fixx){
		printf("Laplacian_Sigma:: Error: x does not match. input=%.3e internal x= %.3e diff = %.3e\n",x,*fixx, x-*fixx);
	}
	double val = 0;

//#if (LAPLACIAN==1||R_FORMULA==1)
	switch(mode){
		case 'l'://regular
#if IBP==1
			val=gsl_spline_eval_deriv(spline_ptr, r,r_accel_ptr);
#if NS==2
			val*=( kt*r*std::cyl_bessel_j(1,r*kt) + 2*pow(r/ns_pow,2)*std::cyl_bessel_j(0,r*kt));
#else
			val*= kt*r*std::cyl_bessel_j(1,r*kt);
#endif

#elif IBP==2
			val=gsl_spline_eval(spline_ptr, r,r_accel_ptr);
			val*=-kt*kt*r*std::cyl_bessel_j(0,r*kt);
#else
			val=gsl_spline_eval_deriv2(spline_ptr, r,r_accel_ptr);
			val+=gsl_spline_eval_deriv(spline_ptr, r,r_accel_ptr)/r;
			val*=r*std::cyl_bessel_j(0,r*kt);
#endif			
			break;
		case 's'://if grid is of laplacian sigma
			val=gsl_spline_eval(spline_ptr, r,r_accel_ptr);
			val*=r*std::cyl_bessel_j(0,r*kt);
			break;
		case 'w'://for weizsacker-william 
			val=gsl_spline_eval(spline_ptr, r,r_accel_ptr);
			val*=std::cyl_bessel_j(0,r*kt)/r;
			break;
		default:
			printf("unknown option in laplacian sigma %c\n",mode);
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
#if NS>=1
	val*=exp(-pow(r/ns_pow,2));
#endif
#if (R_CHANGE_VAR==1)
	return(val/pow(1-rho,2));
#elif (HANKEL==1||R_CHANGE_VAR==0)
	return val;
#endif
}

template <typename Sig>double Laplacian_Sigma<Sig>::constant(double r , const std::vector<double> &par)const {
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
#if NS>=1
	val*=exp(-pow(r/ns_pow,2));
#endif
	return(val);
}


