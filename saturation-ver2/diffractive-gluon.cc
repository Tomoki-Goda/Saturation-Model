#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include<getopt.h>

#include"/home/tomoki/Numerics/clenshaw-curtis.hh"



Clenshaw_Curtis integrator_phi2(64);

double sigma(double x,double r){
	double s0=1/*23*/,x0=0.0003,lam=0.3;
	double Qs2=pow(x0/x,lam);
	double val=s0*(1-exp(-r*r*Qs2/4));
	return(val);
}
/*
double phi2_integrand(double *R,void* v_args){
	double r=*R;
	double* args=(double *)v_args;
	double k=args[0],beta=args[1],xp=args[2];
	double val=r*std::cyl_bessel_j(2,k*r)*std::cyl_bessel_k(2,std::sqrt(beta/(1-beta))*k*r) *sigma.(xp,r);
	return(val);	
}

double phi2(double* vars){
	//double var[]={k,beta,xp};
	double k=*vars;
	double val=integrator_phi2.integrate(&phi2_integrand,(void*)vars, 1.0e-10,100,1.0e-6);
	val*=k*k;
	return(val);
}
*/
double phi2_integrand(double *U,void* v_args){
	double u=*U;
	double* args=(double *)v_args;
	double k=args[0],beta=args[1],xp=args[2];
	double val=u*std::cyl_bessel_j(2,u)*std::cyl_bessel_k(2,std::sqrt(beta/(1-beta))*u) *sigma(xp,u/k);
	return(val);	
}

double phi2(double* vars){
	//double var[]={k,beta,xp};
	double k=*vars;
	double val=integrator_phi2.integrate(&phi2_integrand,(void*)vars, 1.0e-10,100,1.0e-6);
	//val*=k*k;
	return(val);
}

Clenshaw_Curtis integrator_gD(32);

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
	double val=81.0/(256*pow(PI,4)*5.6)*(beta*beta/pow(1-beta,3));
	val*=integrator_gD.integrate(&gD_integrand,(void* )args, 1.0e-10,100,1.0e-5);
	//std::cout<<val<<std::endl;
	return(val*pow(23.0,2)/pow(0.389,2));
}
/////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////

double glu(double k,double Qs2){
	double Qsa2=9.0/4.0*Qs2;
	double val=(23.0/(2*0.389))/(pow(2*PI,2)*PI*Qsa2)*exp(-pow(k,2)/Qsa2);
	return val;
}
Clenshaw_Curtis integrator_k1(32);
double k1_integrand(double*K, void* v_args){
	double k1=*K;
	double *args=(double *)v_args;
	double Qs2=args[2],beta=args[0],k=args[1];
	double val=glu(k1,Qs2);
	val*=pow(beta,2)+pow(1-beta,2)+pow(k1/k,2)*(1-beta)-(pow( (1-2*beta)*pow(k,2) - (1-beta)*pow(k1,2), 2) + 2*beta*(1-beta)*pow(k,4))/(pow(k,2)*pow( pow(pow(k,2)+(1-beta)*pow(k1,2) ,2)-4*pow((1-beta)* k*k1,2),0.5));
	//std::cout<<val<<std::endl;
	return(2*k1*val);
}

Clenshaw_Curtis integrator_k(32);
double dfdY_integrand(double *k,void* v_args ){

	double *args=(double *)v_args;
	double beta=args[0];
	args[1]=*k;

	double val=integrator_k1.integrate(&k1_integrand,(void*)args,1.0e-10,100,1.0e-3 );
	val=val*val;
	//std::cout<<val<<std::endl;
	val*=*k*2*PI;
	return(val);
}

double D_g(double beta){
	double Qs2=0.25;
	double args[]={beta,0,Qs2};
	double val=pow(2*PI,4)/(pow(23.0/0.389,2)*9*Qs2);

	val*=PI/pow(1-beta,3);
	val*=integrator_k.integrate(&dfdY_integrand,(void*)args,1.0e-10,100,1.0e-2);
	return(val);
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////

double gluF(double k,double Qs2){
	double Qsa2=9.0/4.0*Qs2;
	//double Qsa2=Qs2;
	double val=(23.0/0.389)*(4.0/3.0)/(Qsa2*pow(2*PI,2))*exp(-pow(k,2)/Qsa2);
	return(val);
}
double gluF_gbw(double k,double Qs2){
	//double Qsa2=Qs2;
	//double Qsa2=Qs2;
	double val=(23.0/0.389)*(4.0/3.0)/(Qs2*pow(2*PI,2))*exp(-pow(k,2)/Qs2);
	return(val);
}
double k1F_integrand(double*K, void* v_args){
	double k1=*K;
	double *args=(double *)v_args;
	double Qs2=args[2],beta=args[0],k=args[1];
	//double val=gluF(k1,Qs2);
	double val=gluF_gbw(k1,Qs2);

	val*=pow(beta,2)+pow(1-beta,2)+pow(k1/k,2)*(1-beta)-(pow( (1-2*beta)*pow(k,2) - (1-beta)*pow(k1,2), 2) + 2*beta*(1-beta)*pow(k,4))/(pow(k,2)*pow( pow(pow(k,2)+(1-beta)*pow(k1,2) ,2)-4*pow((1-beta)* k*k1,2),0.5));
	//std::cout<<val<<std::endl;
	return(2*k1*val);
}

double dfdYF_integrand(double *k,void* v_args ){
	double *args=(double *)v_args;
	double beta=args[0];
	args[1]=*k;
	double val=integrator_k1.integrate(&k1F_integrand,(void*)args,1.0e-10,100,1.0e-3 );
	val=val*val;
	//std::cout<<val<<std::endl;
	val*=*k*2;
	return(val);
}


double dfdY(double beta,double xp){
	//double Qs2=1*2 ;
	double s0=23.0/0.389,x0=0.0003,lam=0.3;
	double Qs2=pow(x0/xp,lam)/4;

	double args[]={beta,0,Qs2};
	double val=9.0/(64.0*pow(1-beta,3) );

	val*=integrator_k.integrate(&dfdYF_integrand,(void*)args,1.0e-10,100,1.0e-2);
	//val*=(pow(2*PI,4)) / (s0*s0*9*Qs2) ;
	val*=(3*PI)/(2*5.6);
	return(val);
}



////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////



int main(int argc, char** argv){
	static struct option long_options[]={
		{"out",required_argument,0,'o'}
	};
	int opt_ind;
	int c;
	std::string outfile="./points.txt"; 
	while(1){
		c=getopt_long(argc,argv,"o:",long_options,&opt_ind);
		if(c==-1){
			break;
		}
		switch(c){
			case 'o':
				outfile=optarg;
				break;
			default:
				std::cout<<"Error"<<std::endl;
		}
	}

	double p,val;
	std::fstream file;
	file.open(outfile,std::fstream::out);
	for (int i=0;i<5;i++){
		p=0.01+i*(0.99-0.01)/4;
		//val=gD(p,0.0042);
		val=dfdY(p,0.0042);
		//val=D_g(p);
		std::cout<<val<<std::endl;
		file<<p<<"\t"<<val<<std::endl;

	}
	file.close();
	return(0);

}
