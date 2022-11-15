#include<iostream>
#include<cmath>
#include<fstream>

#include"/home/tomoki/Numerics/clenshaw-curtis.hh"



Clenshaw_Curtis integrator_phi2(64);

double sigma(double r,double x){
	double s0=1/*23*/,x0=0.0003,lam=0.3;
	double Qs2=pow(x0/x,lam);
	double val=s0*(1-exp(-r*r*Qs2/4));
	return(val);
}

double phi2_integrand(double *R,void* v_args){
	double r=*R;
	double* args=(double *)v_args;
	double k=args[0],beta=args[1],xp=args[2];
	double val=r*std::cyl_bessel_j(2,k*r)*std::cyl_bessel_k(2,std::sqrt(beta/(1-beta))*k*r) *sigma(r,xp);
	return(val);	
}

double phi2(double* vars){
	//double var[]={k,beta,xp};
	double k=*vars;
	double val=integrator_phi2.integrate(&phi2_integrand,(void*)vars, 1.0e-10,100,1.0e-6);
	val*=k*k;
	return(val);
}

Clenshaw_Curtis integrator_gD(16);

double gD_integrand(double *K, void* v_args){
	double k=*K;
	double *args=(double*)v_args;
	args[0]=k;
	double val=phi2(args);
	val=val*val;
	val*=2 *k;
	return(val);
}
double gD(double beta,double xp){
	double args[3]={0,beta,xp};
	double val=81.0/(256*pow(PI,4)*5.6*xp)*(beta/pow(1-beta,3));
	val*=integrator_gD.integrate(&gD_integrand,(void* )args, 1.0e-10,100,1.0e-4);
	//std::cout<<val<<std::endl;
	return(val/*/pow(0.389,2)*/);
}

int main(){

	double p,val;
	std::fstream file;
	file.open("./points.txt",std::fstream::out);
	for (int i=0;i<10;i++){
		p=0.05+i*(0.9-0.05)/9;
		val=gD(p,0.0042);
		std::cout<<val<<std::endl;
		file<<p<<"\t"<<val<<std::endl;

	}
	file.close();
	return(0);

}
