//static double PSI[5][2*N_SIMPS_R+1][MAXN];//pre-evaluated sets of psi
//static double SAMPLES[5][2*N_SIMPS_R+1][MAXN];//sampled points of integrand 
#define MAXN 600
#define RESCALE 1.05


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
//static double PSI[5][MAXN][2*N_SIMPS_R+1];//pre-evaluated sets of psi
//static double PSI[5][2*N_SIMPS_R+1][MAXN];//pre-evaluated sets of psi
//static double SAMPLES[5][2*N_SIMPS_R +1][MAXN];
static double PSI[5*(2*N_SIMPS_R+1)*MAXN];//pre-evaluated sets of psi
static double SAMPLES[5*(2*N_SIMPS_R +1)*MAXN];
/////////////////////////////////////////////
//static const double ep=R_MIN;//1.0e-6;//for r==0 is divergent or unstable note this value is related to the value chosen for lower limit in chebyshev...
#if (R_CHANGE_VAR==1)
static const double INT_R_MAX=((double)R_MAX)/(1+R_MAX);
static const double INT_R_MIN=((double)R_MIN)/(1+R_MIN);
#else
static const double INT_R_MAX=R_MAX;
static const double INT_R_MIN=R_MIN;
#endif

//static const double R_STEP=r_int_max/(2*N_SIMPS_R);

int N_SIMPS=N_SIMPS_R;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////   generate grid of z-integrated psi values        ////////////////////////////
/////////////////////////////////              for every Q of experimental data   //////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////it writes to global PSI...
void generate_psi_set(){
	double r_step=(INT_R_MAX-INT_R_MIN)/(2*N_SIMPS);
	double r;
	char outline[200];
	sprintf(outline,"r integrated from 0 to %f, with step %.3e. \n\n",INT_R_MAX, r_step);
	log_printf(log_file,outline);
	sprintf(outline,"nf=%d\tN_SIMPS=%d\tN_DATA=%d.\n\n", (int)NF, N_SIMPS,N_DATA);
	log_printf(log_file,outline);
	
	for(unsigned fl=0;fl<(NF-1);fl++){
	
		for(unsigned i=0; i<N_DATA;i++){
			for(unsigned j=0;j<(2*N_SIMPS+1); j++){
				r=r_step*j+INT_R_MIN;
#if (R_CHANGE_VAR==1)
				r=r/(1-r);
				//r=-log(r);
#endif
				//printf("%d\n",fl);
				//*(*(*(PSI+fl )+j )+i )=psisq_z_int(r, *(Q2_DATA+i), fl);
				PSI[i*((NF-1)*(2*N_SIMPS+1)) + j*(NF-1) + fl]=psisq_z_int(r, *(Q2_DATA+i), fl);
			}
		}
	}
#if (PRINT_PROGRESS==1)
	printf("*****************************************\n");
	printf("*            Psi ready                  *\n");
	printf("*****************************************\n");
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// SIMPSON'S INTEGRATION APPROXIMATION /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void sample_integrand(const  double *PSI,  double  *samples, const double* par){
	double xm, r, Q2,term ;
	double r_step=(INT_R_MAX-INT_R_MIN)/(2*N_SIMPS);
	double sudpar[10]={0};
	double sigpar[10]={0};
	parameter(par,sigpar,sudpar);// problem here if -Ofast is used...?
	
	//double ep=1.0e-5;
//////////////////////////////////////////////////////////////
	for(unsigned data_no=0; data_no<N_DATA;data_no++){
		//val=0;
		Q2=Q2_DATA[data_no];

		for(unsigned fl=0;fl<(NF-1);fl++){
			xm=mod_x(X_DATA[data_no], Q2,fl );
			
			for(unsigned j=0;j<(2*N_SIMPS+1); j++){
				
				r=r_step*j+INT_R_MIN;
#if (R_CHANGE_VAR==1)
				r=r/(1-r);
#endif
				//term= (*(*(*(PSI+fl )+j )+data_no )) * ( SIGMA(r,xm,Q2, sigpar,sudpar) )/r;//it should be *r coming from dr r d(theta) but we give r^2 to psi and so /r ;
				term= PSI[data_no*((NF-1)*(2*N_SIMPS+1))+j*(NF-1)+fl ] * ( SIGMA(r,xm,Q2, sigpar,sudpar) )/r;//it should be *r coming from dr r d(theta) but we give r^2 to psi and so /r ;
				//printf("%.2e\n",term);
#if (R_CHANGE_VAR==1)
				term*=pow(1+r,2);
#endif	
				samples[data_no*((NF-1)*(2*N_SIMPS+1))+j*(NF-1)+fl ]=term;
				//samples[fl][j][data_no]=term;
			}
		}
	}

}


void simpson_sum(const double *samples, double * csarray){
	double term,val;
	double r_step=(INT_R_MAX-INT_R_MIN)/(2*N_SIMPS);
	for(int data_no=0;data_no<N_DATA;data_no++){
		val=0;
		for(int j=0;j<(2*N_SIMPS+1);j++){	
			term=0;
			for(int fl=0;fl<(NF-1);fl++){
				//term+=samples[fl][j][data_no];   
				term+=samples[data_no*((NF-1)*(2*N_SIMPS+1))+j*(NF-1)+fl ];
			}	
			if((j==0)||(j==2*N_SIMPS)){
			
			} else if( (j/2)*2==j ){
				term*=2;
			}
			else{
				term*=4;	
			}
			val+=term;
		}
	*(csarray+data_no)=val*(r_step/3);
	}
		
}

int comp_fabs(const void* a, const void* b){
	if(fabs(*(double*)(a))<fabs(*(double*)(b))){
		return(-1);
	}else{
		return(1);
	}
}
double kahn_sum(const double *arr, int len){
//https://en.wikipedia.org/wiki/Kahan_summation_algorithm
	double sum=0 ;
	double c=0;
	
	double y,t;
	for(int i=0;i<len;i++){
		y=arr[i]-c;
		t=sum+y;
		c=(t-sum)-y;
		sum=t;
	}
	return sum;
}
double KBN_sum(const double *arr, int len){
//https://en.wikipedia.org/wiki/Kahan_summation_algorithm
	double sum=0 ;
	double c=0;
	
	double y,t;
	for(int i=0;i<len;i++){
		t=sum+arr[i];
		if(fabs(sum)>=fabs(arr[i])){
			c+=(sum-t)+arr[i];
		}else{
			c+=(arr[i]-t)+sum;
		}
		
		sum=t;
	}
	return (sum+c);
}
void simpson_sum_sorted(const double *samples, double * csarray){
	double summand[(2*N_SIMPS+1)];
	double term,val;

	double r_step=(INT_R_MAX-INT_R_MIN)/(2*N_SIMPS);
	for(int data_no=0;data_no<N_DATA;data_no++){
		val=0;
		for(int j=0;j<(2*N_SIMPS+1);j++){	
			term=0;
			for(int fl=0;fl<(NF-1);fl++){
				//term+=samples[fl][j][data_no];   
				term+=samples[data_no*((NF-1)*(2*N_SIMPS+1))+j*(NF-1)+fl ];
			}	
			if((j==0)||(j==2*N_SIMPS)){
			
			} else if( (j/2)*2==j ){
				term*=2;
			}
			else{
				term*=4;	
			}
			summand[j]=term;
			//val+=term;
		}
		val=0;
		qsort(summand, 2*N_SIMPS+1,sizeof(*summand),&comp_fabs);
		for(int i=0;i<(2*N_SIMPS+1);i++){
			val+=summand[i];
		}
	
		*(csarray+data_no)=val*(r_step/3);
	}		
}

void simpson_kahn_sum(const double *samples, double * csarray){
	double summand[(2*N_SIMPS+1)];
	double term,val;

	double r_step=(INT_R_MAX-INT_R_MIN)/(2*N_SIMPS);
	for(int data_no=0;data_no<N_DATA;data_no++){
		val=0;
		for(int j=0;j<(2*N_SIMPS+1);j++){	
			term=0;
			for(int fl=0;fl<(NF-1);fl++){
				//term+=samples[fl][j][data_no];   
				term+=samples[data_no*((NF-1)*(2*N_SIMPS+1))+j*(NF-1)+fl ];
			}	
			if((j==0)||(j==2*N_SIMPS)){
			
			} else if( (j/2)*2==j ){
				term*=2;
			}
			else{
				term*=4;	
			}
			summand[j]=term;
			//val+=term;
		}
		val=0;
		qsort(summand, 2*N_SIMPS+1,sizeof(*summand),&comp_fabs);
		//val=kahn_sum(summand,2*N_SIMPS+1);
		val=KBN_sum(summand,2*N_SIMPS+1);
	
		*(csarray+data_no)=val*(r_step/3);
	}		
}

void simpson_error(const double *samples, double * error_array){
	double term,val;
	double r_step=(INT_R_MAX-INT_R_MIN)/(2*N_SIMPS);
	double sum[N_DATA*(2*N_SIMPS+1)];//This is large, may hit the limit...
	double fourth[N_DATA];
	for(int data_no=0;data_no<N_DATA;data_no++){
		val=0;
		for(int j=0;j<(2*N_SIMPS+1);j++){	
			term=0;	
			for(int fl=0;fl<(NF-1);fl++){
				term+=samples[data_no*((NF-1)*(2*N_SIMPS+1))+j*(NF-1)+fl ];
			}

			sum[data_no*(2*N_SIMPS+1)+j]=term;	
		}
	}
	FILE* file=fopen("/home/tomoki/Saturation-Model/IntegrandError.txt","w");
	if(file==NULL){
		printf("simpson_error:: file not open\n");
	}//else{
		//printf("FILE OPEN\n");
	//}

	for(int data_no=0;data_no<N_DATA;data_no++){
		fourth[data_no]=0;
		for(int j=0;j<(2*N_SIMPS-3);j++){	
			val= (sum[j]-4.0*sum[j+1]+6.0*sum[j+2]-4.0*sum[j+3]+sum[j+4]);///pow(r_step,4) );
			if(data_no==(N_DATA/2)){
				fprintf(file, "%.5e\t%.5e\t%.5e\n",INT_R_MIN+r_step*(j+2), val,sum[j] );
				//fprintf(file, "%.5e\t%.5e\n",R_MIN+r_step*(j+2),sum[j] );
			};
			val=fabs(val);
			if(val>fourth[data_no]){
				fourth[data_no]=val;
			}
					}
		error_array[data_no]=( fourth[data_no]*(INT_R_MAX-INT_R_MIN) )/180;
		//printf("error : %.3e\n",error_array[data_no]);

	}
	fclose(file);
}
