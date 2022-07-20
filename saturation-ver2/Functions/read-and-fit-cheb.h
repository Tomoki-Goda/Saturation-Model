
#define MAXN 600
#define RESCALE 1.05
//#include"./kahnsum.h"

extern double KBN_sum(const double *arr,int len);
extern double kahn_sum(const double *arr,int len);

double cheb_c(const double * sample_arr, const unsigned* ind1,const unsigned *degree, unsigned dim );
double change_var_revert(double min,double max, double val);
double change_var_revert_log(double min,double max, double val);
void sample(double func(const double *,const double* ), const double * par,  const unsigned * degree, unsigned dim,  double* sample_arr);

extern double SIGMA_PREC;

extern double SIGMA(double , double ,double , const double *, const double*);
extern double psisq_z_int(double, double ,int);
extern double mod_x(double,double, int);

extern void approx_xg(const double *);

extern int parameter(const double*,double*,double*);
///////////Set in main.c//////////////////
extern void log_printf(FILE*,char*);
extern FILE* log_file;

// ////////GLOBAL to this file...////////////////
static double X_DATA[MAXN]={0};
static double Y_DATA[MAXN]={0};
static double wdata[MAXN]={0};
static double Q2_DATA[MAXN]={0};
static double CS_DATA[MAXN]={0};
static double ERR_DATA[MAXN]={0};
static unsigned N_DATA;

//static double FIT_RES[N_PAR+1]={0};

//////////////////GLOBAL ARRAY for DATA/////////////////////

//static double PSI[5*(N_CHEB_R)*MAXN];//pre-evaluated sets of psi
//static double SAMPLES[5*(N_CHEB_R)*MAXN];
/////////////////////////////////////////////
//static const double ep=R_MIN;//1.0e-6;//for r==0 is divergent or unstable note this value is related to the value chosen for lower limit in chebyshev...
//#if (R_CHANGE_VAR==1)
//static const double INT_R_MAX=((double)R_MAX)/(1+R_MAX);
//static const double INT_R_MIN=((double)R_MIN)/(1+R_MIN);
//#else
static const double INT_R_MAX=R_MAX;
static const double INT_R_MIN=R_MIN;
//#endif

extern int N_CHEB;

////////////////////////// FORMAT ////////////////////////////////
double psi_for_cheb(const double *rptr, const double * par){
	double r=change_var_revert(R_MIN,R_MAX,*rptr);
	return psisq_z_int(r, par[0],((int)(par[1]+0.1) ));// *(Q2_DATA+i), fl);
} 

double sigma_for_cheb(const double *rptr,const double *par){
	//printf("%.5e\n",sudpar[0]);
	double r=change_var_revert(R_MIN,R_MAX,*rptr);
	//printf("%.5e\n",par[3]);
 	double val;
 	val= SIGMA(r,par[0] ,par[1], par+2, par+12)/r;
 	
 	//CAUTION!! +12 above is not for good reason, simply par is one array and can't pass array. 
 	//both pars are given 10 places . and 2 for x and Q2. 
 	//AND , /r is r/r^2 for r^2 was given to photon wavefunction to remove singularity.
 	// r comes from jacobian for polar coordinates.
 	return val;
 }
////////////////////////////////////////////////////////////////////


void generate_psi_set(double* psi_arr){
	//double r_step=(INT_R_MAX-INT_R_MIN)/(2*N_SIMPS);
	//double r;
	if(N_CHEB>N_CHEB_R){
		printf("N_CHEB has to be smaller than N_CHEB_R\n");
	}
	char outline[200];
	double par[2];
	double* startptr;
	sprintf(outline,"r integrated from 0 to %f \n\n",INT_R_MAX );
	log_printf(log_file,outline);
	sprintf(outline,"nf=%d\tN_CHEB=%d\tN_DATA=%d.\n\n", (int)NF, N_CHEB ,N_DATA);
	log_printf(log_file,outline);
	
	for(unsigned i=0; i<N_DATA;i++){
		for(unsigned fl=0;fl<(NF-1);fl++){
		
			startptr=psi_arr+(i*(NF-1)*N_CHEB + fl*N_CHEB);
			par[0]=*(Q2_DATA+i);
			par[1]=fl;
			
			sample(&psi_for_cheb,par,&N_CHEB,1,startptr);
		}
	}
#if (PRINT_PROGRESS==1)
	printf("*****************************************\n");
	printf("*            Psi ready                  *\n");
	printf("*****************************************\n");
#endif
}



void sample_integrand(const  double *psi_arr,  double  *samples, const double* par){
	double sigma_arr[N_DATA*(NF-1)*N_CHEB];
	
	double *startptr;
	//double r_step=(INT_R_MAX-INT_R_MIN)/(2*N_SIMPS);
	double sudpar[10]={0};
	double sigpar[10]={0};
	parameter(par,sigpar,sudpar);// problem here if -Ofast is used...?
	double param[22]={0}; 
	///very unfortunately ugly thing to do///
	for(int i=0;i<10;i++){
	
		param[i+2]=sigpar[i];
		param[10+i+2]=sudpar[i];
	//printf(" param: %.3e \t %.3e\n",param[i+2],param[12+i] );
	
	}
	
	//double ep=1.0e-5;
//////////////////////////////////////////////////////////////
	for(unsigned data_no=0; data_no<N_DATA;data_no++){
		//val=0;
		param[1]=Q2_DATA[data_no];

		for(unsigned fl=0;fl<(NF-1);fl++){
		
			startptr=sigma_arr+(data_no*(NF-1)*N_CHEB + fl*N_CHEB);
			
			param[0]=mod_x(X_DATA[data_no], Q2_DATA[data_no]  ,fl );
			sample(&sigma_for_cheb, param, &N_CHEB,1,startptr);		
			
		}
	}
	for(int i=0; i<=(N_DATA*(NF-1)*N_CHEB); i++){
		samples[i]=psi_arr[i]*sigma_arr[i];
	}

}

double curtis_clenshaw_sum(const double *sample_arr){
	double summand[N_CHEB/2+1];
	int ind=0;
	double val;

	
	summand[0]= cheb_c(sample_arr, &ind , &N_CHEB , 1 );
	
	for(int i=1; i<(N_CHEB/2 ) ;i++){
		ind=2*i;
		val=cheb_c(sample_arr, &ind , &N_CHEB , 1 );
		val*=2.0/( 1-ind*ind );
		summand[i]=val;
	}
	val =cheb_c(sample_arr, &ind , &N_CHEB , 1 ) ;
	val*=1.0/(1-N_CHEB*N_CHEB);
	summand[N_CHEB/2]=val;
	
	
	return KBN_sum(summand,N_CHEB/2+1);
	 
}

void generate_data_set(const double *par,const double *psi_arr, double *csarray){
	if(2*(N_CHEB/2)!=N_CHEB){
		printf("N_CHEB has to be multiple of 2\n");
	}
	//double psi_arr[(NF-1)*N_CHEB*N_DATA];
	static double samples[5*N_CHEB_R*MAXN];
	
	sample_integrand(psi_arr,samples ,par);
	
	double val;
	int shift;
	for(unsigned data_no=0; data_no<N_DATA;data_no++){
		val=0;
		for(unsigned fl=0;fl<(NF-1);fl++){
			shift=(data_no*(NF-1)*N_CHEB + fl*N_CHEB);
			 val+=curtis_clenshaw_sum(samples+shift );
			
		}
		csarray[data_no]=val*(R_MAX-R_MIN)/2;
	}
	//simpson_sum(SAMPLES, csarray);	
	//simpson_sum_sorted(SAMPLES, csarray);//intention of this is to add small values first to avoid loss by rounding.	
	
//////////////////////////////////////////////////////////////
}




















