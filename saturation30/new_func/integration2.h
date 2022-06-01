//#include<stdio.h>
//#include<math.h>

//#define FLAVOUR 2 
#if MODEL==0
#define SIGMA sigma_gbw
#elif MODEL==1
#define SIGMA sigma_bgk
#elif MODEL==2
#define SIGMA sigma_gbs
#endif

#define MAXN 600
#define N_SIMPS_R 150
//////////////////GLOBAL ARRAY for DATA/////////////////////
static double PSI[NF-1][MAXN][2*N_SIMPS_R+1];//pre-evaluated sets of psi
static double *X_VALS;//[N_DATA];
static double *Q2_VALS;//[N_DATA];
/////////////////////////////////////////////
static const double r_int_max=30.0;
static const double r_step=r_int_max/(2*N_SIMPS_R);
static const double ep=0.0;//for r==0 is divergent or unstable


void import_points(double * X_DATA,double *Q2_DATA){
	//because they are static in other file...
	X_VALS=X_DATA;
	Q2_VALS=Q2_DATA;
	
}

double mod_x(double x, double Q2, unsigned flavour) {
	double m_fsq;

	switch (flavour) {
	case 0:
		m_fsq = MASS_L2;
		break;
	case 1:
		m_fsq = MASS_S2;
		break;
	case 2:
		m_fsq = MASS_C2;
		break;
	case 3:
		m_fsq = MASS_B2;
		break;
	default:
		printf("wrong input %c\n",flavour);
		m_fsq = MASS_L2;
	}
	return (x * (1.0 +( 4.0 * (m_fsq/Q2)) ));
}

////////////////////////////////generate grid of z-integrated psi values//////////////////////////////////
/////////////////it writes to global PSI...
void generate_psi_set(){
	for(unsigned fl=0;fl<(NF-1);fl++){
		for(unsigned i=0; i<N_DATA;i++){
			for(unsigned j=0;j<(2*N_SIMPS_R+1); j++){			
				*(*(*(PSI+fl )+j )+i )=psisq_z_int_double(R_STEP*j+ep, *(Q2_VALS+i), fl+0.5);
			}
		}
	}
}

/////////////////////////////now integrate over r////////////////////////////

void generate_data_set(double *par double *csarray){
	//csarray is counterpart of CS_DATA ...
	//double integral[N_DATA];
	double term;
	double r,Q2,xm;
	for(unsigned i=0; i<N_DATA;i++){
		val=0;
		for(unsigned fl=0;fl<(NF-1);fl++){
			for(unsigned j=0;j<(2*N_SIMPS_R+1); j++){	
				r=R_STEP*j+ep;
				Q2=*(Q2_VALS+i);
				xm=mod_x(*(X_VALS+i), Q2,fl );
				
				term= (*(*(PSI+fl )+j )) * ( SIGMA(r,xm, par) )/r//it should be *r coming from dr r d(theta) but we give r^2 to psi and so /r ;
				if((j==0)||(j==2*N_SIMPS_R)){
					
				} else if( (j/2)*2==j ){
					term*=2;
				}
				else{
					term*=4;	
				}
				
			}
		}
		*(csarray+i)=term*(R_STEP/3);
		
	}
}

