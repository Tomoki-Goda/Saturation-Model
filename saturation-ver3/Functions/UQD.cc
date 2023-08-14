#include<cmath>
#include<iostream>
#include<vector>
#include<string>
#include<complex>
#include<chrono>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_spline.h>
#include<gsl/gsl_errno.h> 
#include<gsl/gsl_interp2d.h>
#include<gsl/gsl_spline2d.h>
#include<cuba.h>
//#include"miscellaneous.hh"

static double change_var(double & var,double &  jac,const double min, const double max,const double c){//This version is (in theory) regular at max->Inf
	double den=( (c==1)?(1):(c+var*(1-c)) );
	jac= ( (min==0.0)?(c*pow(den,-2)*max):(c*pow(den,-2)*(max-min)) ) ;
	var= ( (min==0.0)?(max*var):((max*var+c*min*(1-var)) ))/den;
	//var= (max*var+min*c*(1-var))/den;
	
//#if TEST==1	
	if(var>max) {
		if(fabs((var-max)/max)>1.0e-15){
			printf("value below limit %.3e -> %.3e [%.3e, %.3e] diff %.3e, c=%.3e\n",(1-den)/(1-c),var,min,max,var-max, c);
		}
		var=max;
	}else if(var<min){
		if(fabs((min-var)/min)>1.0e-15){
		printf("value below limit %.3e -> %.3e [%.3e, %.3e] diff %.3e, c=%.3e\n",(1-den)/(1-c),var,min,max,min-var, c);
		}
		var=min;
	}
//#endif
	return var;
}


class unintegrated_quark{
	private:
	gsl_interp_accel *x_accel_ptr, *k2_accel_ptr;
	gsl_spline2d *  spline_ptr;
	double x, Q2, d2, k2_min, k2_max;
	double as=0.2,nf=0.5,PI=3.14159265359;
	//nf=0.5 is because we want q(x,k,mu) 
	//see Eellis, Webber, Sterling Eq.4.99
	inline double split_qg(double z, double kd);
	
	public:
	int k2_npts,x_npts;
	double k2_array[250],x_array[250],aF_array[62500];
	double integrand( double z, double k2);
	void load(std::string s);
	inline double gluon(double x, double k2);
	inline double  operator()( double Q2,double z, double k2);
	unintegrated_quark(){
		x_accel_ptr = gsl_interp_accel_alloc ();
		k2_accel_ptr = gsl_interp_accel_alloc ();	
	}
	~unintegrated_quark(){
		gsl_spline2d_free (spline_ptr);
		gsl_interp_accel_free (x_accel_ptr);
		gsl_interp_accel_free (k2_accel_ptr);	
	}

};


void unintegrated_quark::load(std::string s){
//	double k2_array[250],x_array[250],aF_array[62500];
	FILE* file=fopen(s.c_str(),"r");
	printf("%s\n",s.c_str());
	
	double x,k2,f;
	double x_p,k2_p;
	int xi=0, k2i=0,fi=0,kflag=0;
	for(int i=0;i<62500;i++){
		fscanf(file,"%lf %lf %lf",&x,&k2,&f);
		if(feof(file) ||(i+1==62500)){
			break;
		}
		
		
		x=exp(x);
		k2=exp(k2);
		aF_array[fi++]=f;
/////////Record x ////////////////////		
		if(i==0){
			x_p=x;
			x_array[xi++]=x;
		}else{
			if( fabs(x_p-x)/fabs(x_p+x)>1.0e-10 ){//only whan x changes
				//printf("x=%.3e x_prev=%.3e\n",x,x_p);
				x_p=x;
				x_array[xi++]=x;
				kflag=-1;//stop recording k2
			}
		}
		
////////Record k2 //////////////////
		if(kflag!=-1){
			k2_array[k2i++]=k2;
		}		
	}
	fclose(file);
	
	k2_npts=k2i;
	x_npts=xi;
	
	k2_min=k2_array[0];
	k2_max=k2_array[k2_npts-1];
	printf("k2: %d [%.2e, %.2e], x: %d [%.2e, %.2e]\n",k2_npts,k2_array[0],k2_array[k2_npts-1], x_npts,x_array[0],x_array[x_npts-1]);
	//for(int i=0;i<x_npts;i++){
	//	printf("%.2e \n" ,x_array[i]);
	//}printf("\n");
	//for(int i=0;i<k2_npts;i++){
	//	printf("%.2e \n" ,k2_array[i]);
	//}
	spline_ptr = gsl_spline2d_alloc(gsl_interp2d_bicubic,k2_npts, x_npts);
	gsl_spline2d_init (spline_ptr,k2_array, x_array, aF_array, k2_npts, x_npts);	
}
inline double unintegrated_quark::gluon(double x, double k2){
	return	gsl_spline2d_eval_extrap(spline_ptr ,k2, x, k2_accel_ptr, x_accel_ptr);
 	//if(x<x_array[0]||x>x_array[x_npts-1]||k2<k2_array[0]||k2>k2_array[k2_npts-1]){
 	//	printf("taking x= %.2e [%.2e, %.2e], k2=%.2e ,[%.2e, %.2e] \n",x,x_array[0],x_array[x_npts-1],k2,k2_array[0],k2_array[k2_npts-1]);
 	//	return 0;
 	//}
 	//return	gsl_spline2d_eval(spline_ptr ,k2, x, k2_accel_ptr, x_accel_ptr);
 	
}

inline double unintegrated_quark::split_qg(double z, double kd2){
	//splitting function Pqg in hep-ph/9405388
	// there is typo in 2110.06156
	double result;
	double z2=pow(z,2),z12=pow(1-z,2);
	result=(as*nf)/(2*PI)*pow((1+z*(1-z)*kd2),-2)*(z2+z12+4*z2*z12*kd2);
	return(result);
}
inline double unintegrated_quark::integrand( double x1, double x2){
	double result, jac1,jac2, k2=x2,z=x1;
	change_var(k2,jac1,k2_min, k2_max,fabs(k2_max-k2_min));
	change_var(z,jac2,this->x, 1 ,1-(this->x));
	//printf("k2: %.2e->%.2e  [%.2e, %.2e]\n",x2,k2,k2_min, k2_max);
	result= ((Q2-d2/(1-z)-z*k2)<0)?(0):(jac1*jac2*split_qg(z,k2/d2)*gluon(x/z,k2));
	//result= ((Q2-d2/(1-z)-z*k2)<0)?(0):(jac1*jac2*split_qg(z,k2/d2)*gluon(x,k2));
	if(std::isnan(result)){
		printf("res=%.2e, z=%.2e k2=%.2e, Q2=%.2e x=%.2e Delta^2=%.2e\n",result,z,k2,Q2,x,d2);
	}
	return result;
}
int cuba_integrand(const int * ndim, const double  * intv,const int * ncomp,double *  f, void* p){
	unintegrated_quark* intg= (unintegrated_quark*)(p);
	*f=intg->integrand(intv[0],  intv[1]);
	return 0;
}
inline double  unintegrated_quark::operator()( double x, double d2, double Q2){
			const long long int mineval=1000, maxeval=1000000;//use llChure if larger than ~1.0e+9
			const long long int nstart=100,nincrease=100;
			long long int neval=0;
			const int flag= 0+4*0+8*1+16*0+32*0;
			int nregions=0,fail=0;
			double  integral[3]={0},error[3]={0},prob[3]={0};
			char statefile[100]="";
			double  result=0;
			int ndim=2,key=13;
			
			this->Q2=Q2;
			this->x=x;
			this->d2=d2;
			
			llCuhre(ndim, 1,
				&cuba_integrand,
				(void*)this,
				 1,1.0e-3 ,1.0e-5, flag, mineval,maxeval, key,NULL,NULL, &nregions, &neval,  &fail, integral, error, prob
			);
			
			return(integral[0]);
}

int main(int argc, char** argv){
 unintegrated_quark sigma;
 sigma.load(argv[1]);
 double x=0 ,k2=0 ,q2=100;
 FILE* file=fopen(argv[2],"w");
 


	for(int i=0;i<sigma.x_npts;i++){
		x=sigma.x_array[i];
		for(int j=0;j<sigma.k2_npts;j++){
			k2=sigma.k2_array[j];
		for(int k=0;k<25;k++){
			q2=0.2*pow(1e+5/0.2,((double)k)/2.4e+1);
			fprintf(file,"%.5e\t%.5e\t%.5e\t%.5e\n",log(x),log(k2),log(q2),sigma(x,k2,q2)/k2 );
		}//fprintf(file,"\n" );
	}
 }
 fclose(file);
}



