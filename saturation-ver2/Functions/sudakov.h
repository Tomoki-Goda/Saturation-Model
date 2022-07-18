//#include<stdio.h>
//#include<math.h>


//#include"./control_tmp.h"
//#include"control.h"
//#include"./control-default.h"
//#include"constants.h"

//#include"dipole-cross-section.h"

#include"./mu2.h"


double SIGMA_PREC;
//extern int N_DATA;


//extern double sudakov(double ,double, double* );// in dipole-cross-section
//extern double BASE_SIGMA(double, double, double, double*);
extern double sigma_gbw(double, double, double, const double*);
extern double sigma_bgk(double, double, double, const double*);

///////////////////////////////////////////////////////////////////////////////////////////////
extern double dgauss_(double(* )(double*),double*,double*,double*); 
extern double dgquad_(double(* )(double*),double*,double*,int*);
extern double dadapt_(double(* )(double*),double*,double*,int*,double*,double* ,double*,double*);

static double VAR[3];
//static double DUMMY_ARRAY[20];
static const double *SIGPAR;//=DUMMY_ARRAY;
static const double *SUDPAR;//=DUMMY_ARRAY;

///////////////////////////////////////////////////////////////////////////////////
////////////////////////////    SUDAKOV    ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

double sudakov_p(double mub2, double q2,const double* par) {
	if (q2 <= mub2) {//ensures that lower limit of integral is smaller than upper limit...
		return(0.0);
	}
	if (mub2 < LQCD2){
		printf("\nsudakov:: mu_b is too low!!!\t mub2= %.3e\n",mub2);//\n%f\t%f\n\n",C,r_max);
		getchar();
		
	}
	if (q2 < LQCD2 ) {
		printf("\nsudakov:: Q2 is too low!!!\t Q2= %.3e\t mub2= %.3e\n",q2,mub2);
		//printf("%.3e\t%.3e\t%.3e\n",par[0],pow(r,-2.0),pow(par[1],-2.0));
		getchar();
	}
	double b0 = ((double)(11*CA-2*NF))/12;
	double L_Q_l=log(q2/LQCD2);
	double L_mub_l=log(mub2/LQCD2);
	
	double val=((double)CA)/(2*b0*PI)*(L_Q_l*(log(L_Q_l /L_mub_l)-1) + L_mub_l );
	
	if(val<0){
		if(val<-1.0e-12){
			printf("sudakov p too small: %.3e mub2=%.3e Q2=%.3e\n", val, mub2,q2);
			printf("logQ2/l2 = %.5e logmu2/l2 = %.5e \n", L_Q_l, L_mub_l);
		}
		//getchar();
		val=0;
	}
	if(isnan(val)!=0){
		printf("sudakov p nan: %.3e\n", val);
		getchar();
		
	}
	
	
	return val;
}
double sudakov_np(double r,double mub2, double q2,const double* par) {
	double C= (*(par));
	//double r_max=( *(par+1));
	double g1=(*(par+2));
	double g2=(*(par+3));
	
	if (mub2 < LQCD2){
		printf("\nsudakov:: mu_b is too low!!!\t mub2= %.3e\n",mub2);
	}
	
	double val=(g1/2) *r*r;
	if(q2>(Q0*Q0)){
		val+=(g2/4) * log(mub2*(r*r)/C )*log( q2 /pow(Q0 ,2) );
	}
	
//	if(val<0){
//		if(val<-1.0e-12){
//			printf("sudakov np too small: %.3e r=%.3e  mub2=%.3e Q2=%.3e\n", val,r, mub2,q2);
//		}
//	}
	if(isnan(val)!=0){
		printf("sudakov np nan: %.3e\n", val);
		getchar();
	}
	
	return val;
}

double exp_sud(double r, double mub2, double q2){
	double sud=0.0;
	
	
#if (SUDAKOV>=1)
	sud+=sudakov_p(mub2 ,q2, SUDPAR);
#endif
#if (SUDAKOV==2)
	sud+=sudakov_np(r, mub2 ,q2, SUDPAR);
#endif
	return(exp(-sud));
}
	
///////////////////////////////////////////////////////////////////////////////////
////////////////////////////    JACOBIAN    ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

double dsnpdr(double r, const double * mub2_arr, double q2,const double* sudpar){
#if (SUDAKOV<=1) 
	return 0.0;
#endif 
	//double rmax=sudpar[1];
	double g1=sudpar[2];
	double g2=sudpar[3];
	double val;
	
	double mu2=mub2_arr[0];
	double jacfirst=mub2_arr[1];//rmu2_jac_first(r,sudpar);
	
	val=g1*r;
	if(q2>(Q0*Q0)){
		val += (g2/4)* (2.0/r+jacfirst/mu2)*log(q2/(Q0*Q0));
	}
	
	return(val);
}


double ddsnpdrdr(double r,const double * mub2_arr, double q2,const double* sudpar){
#if (SUDAKOV<=1) 
	return 0.0;
#endif 	
//	double rmax=sudpar[1];
	double g1=sudpar[2];
	double g2=sudpar[3];
	double val;
	double mu2=mub2_arr[0];
	double jacfirst=mub2_arr[1];
	double jacsecond=mub2_arr[2];
	
	val=g1;
	if(q2>(Q0*Q0)){
		val+= (g2/4)*(- 2.0/(r*r) +jacsecond/mu2 -(jacfirst*jacfirst)/(mu2*mu2))*log(q2/(Q0*Q0)) ; 
	}

	return(val);
}

/////////////////////////////////////////////////////////////////////////////////////////

double dspertdr(const double * mub2_arr,double q2, const  double* sudpar ){
#if (SUDAKOV==0)
	return 0.0;
#endif
	double mu2=mub2_arr[0];
	
	if(mu2>=q2){ 
		return(0.0);
	}
	if(mu2<LQCD2){
		printf( "drds:: mu2 out of range");
	}
	double jacfirst=mub2_arr[1];
	double b0= ((double)(11.0*CA-2.0*NF))/12.0;
	double alpha =1.0/(b0* log( mu2/LQCD2) );
	
	double val= - jacfirst* alpha*log(q2/mu2)/mu2;

	val*=((double)CA)/(2.0*PI);

	return(val);
}

double ddspertdrdr(const double * mub2_arr,double q2,const double * sudpar){
#if (SUDAKOV==0)
	return 0.0;
#endif
	double mu2=mub2_arr[0];
		
	if(mu2>q2){
		return(0.0);
	}
	if(mu2<LQCD2){
                printf( "drds:: mu2 out of range");
        }
        double jacfirst=mub2_arr[1];
	double jacsecond=mub2_arr[2];
	
	double b0 = ((double)(11.0*CA-2.0*NF))/12.0;
	double val;

	double logQmu=log(q2/mu2);
	double alpha =1.0/(b0* log( mu2/LQCD2) );

	//derivative s wrt mu2 with jacobians
	double val1	= - jacsecond * alpha*logQmu/mu2;	
	double val2	=   jacfirst*jacfirst * (alpha/pow(mu2,2)) * ( 1 + logQmu*(1+b0*alpha));
	
	val=val1+val2;
	val*= ((double)CA)/(2.0*PI) ; 

	return(val);

}
//////////////////////////////////////////////////////////////////////////
double dsdr(double r, const double * mu2_arr, double q2, const double *SUDPAR){
	double val;
	if((mu2_arr[0]) < q2){
		val=dspertdr(mu2_arr, q2, SUDPAR);
#if (SUDAKOV==2)
		val+=dsnpdr(r,mu2_arr, q2, SUDPAR);
#endif
	}else{
#if (SUDAKOV==2)
		val=dsnpdr(r,mu2_arr, q2, SUDPAR);		
#else
		return(0.0);
#endif
	}
	return(val);
}
double ddsdrdr(double r ,const double * mu2_arr, double q2, const double *SUDPAR){
	double val;
	if((mu2_arr[0]) < q2){
		val=ddspertdrdr(mu2_arr,q2, SUDPAR);
#if (SUDAKOV==2)
		val+=ddsnpdrdr(r,mu2_arr,q2, SUDPAR);
#endif
	}else{
#if (SUDAKOV==2)
		val=ddsnpdrdr(r,mu2_arr,q2, SUDPAR);
#else
		return(0.0);
#endif
	}
	return(val);
}

///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////     Integration     //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
double integrand( double * r_ptr){
	
	double val;
	double r=*r_ptr;
	double  rho=VAR[0];
	double  x=VAR[1];
	double  q2=VAR[2];
	
	double dels,deldels;
	
	double mu2_arr[3]={0};

	int signal=compute_mu2(r, SUDPAR, mu2_arr,3);//compute mu2, dmu2/dr dd mu2/drdr
///////////return 0 for invalid par/////////////
	if(signal!=0){
		printf("integrand:: mu2= %.3e %.3e %.3e r=%.3e  C=%.3e rmax=%.3e\n",mu2_arr[0] ,mu2_arr[1] ,mu2_arr[2] ,r,SUDPAR[0],SUDPAR[1] );
		return 0;
	}
#if ((MODEL==22)||(MODEL==2))
	if(SIGPAR[2]<0){
		return 0;//also delt with in the sigma_gbw
	}
#endif
#if ((MODEL==3)||(MODEL==1))
	if(SIGPAR[3]<0){
		return 0;
	}
#endif
////////////////////////////////////////////////
	
	if((mu2_arr[0]) <LQCD2){
		printf("integrand::mu2 out of range\n");
		printf("R %f  Q2 %f x %f r %f mu2 %f\n" , rho, q2, x, r , (mu2_arr[0]) );
		getchar();
	}

#if (SUDAKOV==1)	
	if((mu2_arr[0]) > q2){
		printf("integrand:: mu2=%.3e Q2=%.3e\n",mu2_arr[0], q2 );
		return(0.0);
	}
#endif
	dels=dsdr(r, mu2_arr, q2, SUDPAR);
	deldels=ddsdrdr(r,mu2_arr, q2, SUDPAR);
	
	
	double logrr=log(rho/r);

	val=(2.0-logrr)*dels+r*logrr*(dels*dels-deldels) ;
	val=val * exp_sud( r , mu2_arr[0], q2 ) * BASE_SIGMA(r,x,q2,SIGPAR);
	
	if(isnan(val)!=0){
		printf("\n***************INTEGRAND*********************\n");
		printf("value= %f for r= %.2e \n",val,r);
		printf("rho=%.2e Q2=%.2e x=%.2e\n",rho ,q2,x);
		printf("ds/dr=%.2e dds/drdr= %.2e log(R/r) = %.2e\n",dels,deldels,logrr);
		printf("sudakov= %.2e sigma=%.2e\n", exp_sud( r , mu2_arr[0], q2 ) ,BASE_SIGMA(r,x,q2,SIGPAR));
		
		printf("Parameters\n");
		printf("SIGPAR: %f %f %f %f\n",SIGPAR[0],SIGPAR[1],SIGPAR[2],SIGPAR[3]);
#if (SUDAKOV==1)  
		printf("SUDPAR: %f %f\n",SUDPAR[0],SUDPAR[1]);
#elif (SUDAKOV==2)
		printf("SUDPAR: %f %f %f %f \n",SUDPAR[0],SUDPAR[1],SUDPAR[2],SUDPAR[3] );
#endif
		printf("mu2 dmu2dr ddmu2/drdr %f %f %f \n", mu2_arr[0],mu2_arr[1],mu2_arr[2]);
		//getchar();
		return 0;
	}

	return val;
}




double integral_term(double r, double x, double q2,const  double * sigmapar,const  double* sudpar){
	double result=0;
	double rmin=R_MIN;
	
#if (SUDAKOV<=1)
	double rmin_2;
	int signal=rmin2(q2, SUDPAR,&rmin_2 );
	if(signal==0){
		if(rmin<0){
			printf("error rmin neg\n");
		}
		rmin=sqrt(rmin_2);
		if(rmin>r){
			return 0.0;
		}
	}else if(signal==9){
		return 0.0;
	}else if(signal==1){
		printf("C is negative\n");
		return 0.0;
	}else{
		printf( "unrecognized signal from rmin2\n");
	}
#endif

	////////////////////////////////////////////
	//int N=96;
	//result=dgquad_(&integrand,&rmin,VAR,&N);
	///////////////////////////////////////////
	//double N=DGAUSS_PREC*0.01;
	//result=dgauss_(&integrand,&rmin,VAR,&N);
	////////////////////////////////////
	int seg=1;
	
	double NRel=SIGMA_PREC ; //SIGMA_PREC is global and is controled in main.c or read-and-fit.c 
	
	//double NRel=DGAUSS_PREC ;
	double NAbs=1.0e-10;
	double error=0;
	//printf("rmin= %f\trmax =%f\n",rmin, *VAR);
	dadapt_(&integrand,&rmin,VAR,&seg ,&NRel, &NAbs, &result, &error)	;
	
	
////////////////////////////////////////////////////////////////////////////////
	if(isnan(result)!=0){
		printf("***************INTEGRAL*********************\n");
		printf("%f \t %f %f\n",result, rmin,VAR[0]);
		printf("Parameters\n");
		printf("%f %f %f %f\n",SIGPAR[0],SIGPAR[1],SIGPAR[2],SIGPAR[3]);
#if (SUDAKOV==1)  
		printf("%f %f\n",SUDPAR[0],SUDPAR[1]);
#elif (SUDAKOV==2)
		printf("%f %f %f %f \n",SUDPAR[0],SUDPAR[1],SUDPAR[2],SUDPAR[3] );
#endif
		//printf("%f %f %f \n", mu2_arr[0],mu2_arr[1],mu2_arr[2]);
		//getchar();
		return 0;
	}
////////////////////////////////////////////////////////////////////////////////////

	return result;
	
}

double sigma_s(double r, double x, double q2, const double * sigmapar, const double* sudpar){
	
	VAR[0]=r;
	VAR[1]=x;
	VAR[2]=q2;
	SIGPAR=sigmapar;
	SUDPAR=sudpar;
///////////return 0 for invalid par/////////////
	//if((sudpar[0]<0)||(sudpar[1]<0)){
	if(sudpar[0]<0){
		return 0;
	}
	if(sudpar[1]<Q0 ){
		return 0;
	}
	
#if (SUDAKOV==2)
	if((sudpar[2]<0)||(sudpar[3]<0)){
		return 0;
	}
#endif
	double val=BASE_SIGMA(r,x,q2, sigmapar );
#if (SUDAKOV==0)
	return val;
#endif
//////////////////////////////////////////////
	
	double sud=0.0;
	double mu2;
	int signal=compute_mu2(r, SUDPAR, &mu2,1);//compute mu2, dmu2/dr dd mu2/drdr
///////////return 0 for invalid par/////////////
	if(signal!=0){
		printf("sigma_s:: mu2= %.3e r=%.3e  C=%.3e rmax=%.3e\n",mu2,r,SUDPAR[0],SUDPAR[1] );
		getchar();
		return 0;
	}
////////////////////////////////////////////////
	
#if (SUDAKOV==1)  	
	if(mu2>q2){
		return val;
	}
#endif	
	val*=exp_sud(r,mu2,q2);

	val+=integral_term(r,x,q2,sigmapar,sudpar);
	
////////////////////////////////////////////////////////////////////////////////////	
	if(isnan(val)!=0){
		printf("%f %f\n",val,sud);
		printf("Parameters\n");
		printf("%f %f %f %f\n",sigmapar[0],sigmapar[1],sigmapar[2],sigmapar[3]);
#if (SUDAKOV==1)  
		printf("%f %f\n",sudpar[0],sudpar[1]);
#elif (SUDAKOV==2)
		printf("%f %f %f %f \n",sudpar[0],sudpar[1],sudpar[2],sudpar[3] );
#endif
		printf("mu2= %f\n",mu2);
		return 0;
		//getchar();
	}
////////////////////////////////////////////////////////////////////////////////////
	return val;
}





