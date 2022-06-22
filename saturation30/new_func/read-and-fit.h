#define MAXN 600
#define RESCALE 1.05



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

static double FIT_RES[N_PAR];

//////////////////GLOBAL ARRAY for DATA/////////////////////
//static double PSI[5][MAXN][2*N_SIMPS_R+1];//pre-evaluated sets of psi
static double PSI[5][2*N_SIMPS_R+1][MAXN];//pre-evaluated sets of psi
/////////////////////////////////////////////
static const double ep=1.0e-5;//for r==0 is divergent or unstable note this value is related to the value chosen for lower limit in chebyshev...
#if R_CHANGE_VAR==1
static const double r_int_max=0.97;
#else
static const double r_int_max=30.0;
#endif

static const double R_STEP=r_int_max/(2*N_SIMPS_R);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////   generate grid of z-integrated psi values        ////////////////////////////
/////////////////////////////////              for every Q of experimental data   //////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////it writes to global PSI...
void generate_psi_set(){
	double r;
	char outline[200];
	sprintf(outline,"r integrated from 0 to %f, with step %.3e. \n\n", r_int_max, R_STEP);
	log_printf(log_file,outline);
	sprintf(outline,"nf=%d\tN_SIMPS=%d\tN_DATA=%d.\n\n", (int)NF, N_SIMPS_R,N_DATA);
	log_printf(log_file,outline);
	
	for(unsigned fl=0;fl<(NF-1);fl++){
		for(unsigned i=0; i<N_DATA;i++){
			for(unsigned j=0;j<(2*N_SIMPS_R+1); j++){
				r=R_STEP*j+ep;
#if R_CHANGE_VAR==1
				r=r/(1-r);
#endif
				*(*(*(PSI+fl )+j )+i )=psisq_z_int(r, *(Q2_DATA+i), fl);
			}
		}
	}
	printf("*****************************************\n");
	printf("*            Psi ready                  *\n");
	printf("*****************************************\n");
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////  now integrate over r with Simpsons method    /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

void generate_data_set(double *par, double *csarray){
	//csarray is counterpart of CS_DATA ...
	//double integral[N_DATA];
	double term ,val;
	double r,Q2,xm;
	
///////////////////////////unfortunately positions of parameters are now incompatible between models ....///////////////////
	double* sigpar= par;
#if MODEL==22||MODEL==2
	double *sudpar;
	sudpar=(par+3);
#elif MODEL==3
	double sudpar[10];
	sudpar[0]=par[3]*par[5] ;//C*C2
	sudpar[1]=par[4]/par[5];//rmax mu02=C/rmax^2
#if SUDAKOV==2
	sudpar[3]=par[6];
	sudpar[4]=par[7];
#endif
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	for(unsigned i=0; i<N_DATA;i++){
		val=0;
		for(unsigned fl=0;fl<(NF-1);fl++){
			for(unsigned j=0;j<(2*N_SIMPS_R+1); j++){	
				r=R_STEP*j+ep;
#if R_CHANGE_VAR==1
				r=r/(1-r);
#endif
				Q2=*(Q2_DATA+i);
				xm=mod_x(*(X_DATA+i), Q2,fl );

				term= (*(*(*(PSI+fl )+j )+i )) * ( SIGMA(r,xm,Q2, sigpar,sudpar) )/r;//it should be *r coming from dr r d(theta) but we give r^2 to psi and so /r ;
#if R_CHANGE_VAR==1
				term*=pow(1+r,2);
#endif
				//printf("%f\n",term);
				if((j==0)||(j==2*N_SIMPS_R)){
					
				} else if( (j/2)*2==j ){
					term*=2;
				}
				else{
					term*=4;	
				}
				//val+=pow(1-r,-2)*term;
				val+=term;
				
			}
		}
		//printf("%f\n",val);
		*(csarray+i)=val*(R_STEP/3);
		
	}
	//printf("crosssection ready\n");
}


//extern double sigma_DIS( double ,double, double ,double* );
/*******************************************************************************
* The function which returns maximal of two arguments
*******************************************************************************/
double maximum (double arg1, double arg2) {

    double value;
    if(arg1>arg2)
	value = arg1;
    else if(arg1<=arg2)
	value = arg2;
    return value;
}
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////   fitting   ////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void dum_func(void){
	//do nothing
}

double compute_chisq(double *par){
	
	//////////////////////////////////////////
	////////////   GLOBALS   /////////////////
	//////////////////////////////////////////
	//double * CS_DATA	experimental data array
	//double * ERR_DATA	experimental error
	//unsigned N_DATA	number of data 
	/////////////////////////////////////////
	char outline[200];//entry in the log file 
	double chisq=0.0;
	double val;
	clock_t time;
	time=clock();
	double cs[N_DATA];//computed cross-section
	generate_data_set(par,cs);
	
	for(unsigned i=0;i<N_DATA;i++){
		chisq+=pow( ((*(cs+i))-(*(CS_DATA+i)))/(*(ERR_DATA+i)),2);
	}
	time-=clock();
	
	for(unsigned i=0;i<N_PAR;i++){
		sprintf(outline, "%.3e, ",*(par+i));
		log_printf(log_file,outline);
		*(FIT_RES+i)=*(par+i);
	}
	*(FIT_RES+N_PAR)=chisq;
	
	
	sprintf(outline,"%.2e / %d = %.3e, in %.2e sec\n",chisq, N_DATA, chisq/N_DATA, -((double)time)/CLOCKS_PER_SEC);
	log_printf(log_file,outline);
	
	
	return(chisq);
}

void fcn(int npar, double grad[], double*fcnval, double *par,unsigned iflag,void (*dum)(void) ){
	/////////////////////////////////////////
	//for detail see MINUIT documentatin.
#if (MODEL==1||MODEL==3)	
	approx_xg(par+1);//generate chebyshev coefficients
#endif
	*fcnval=compute_chisq(par);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////        READ           //////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int load_data(){
	double dum, sysup, sysdown, stat, fac;
	
	unsigned i=0;
	unsigned j=0;
	FILE* file;
	/////////////////////////////// HERA /////////////////////////////////////
	//fprintf(stdout, "HERA tot\n");
	fprintf(stdout, "HERA tot\n");
	file=fopen("../data/hera_tot.dat","r");
	
	double alpha =1.0/137 ;//fine structure const 1/137;
	double xmp0 = 0.93827;//proton mass in GeV
	double units =1.0/389.40; //2.56819e-3; //micro-barn to GeV^-2
	j=0;
	while((!feof(file))&&(j<597)){
		fscanf(file,"%lE %lE %lE %lE %lE", (Q2_DATA+i),(X_DATA+i),(wdata+i),(CS_DATA+i),(ERR_DATA+i)); 
		/////formula in I. Abt et al 2017////////
		fac = pow(Q2_DATA[i],2)*(1-X_DATA[i])/ (4*pow(PI,2)*alpha*(Q2_DATA[i]+pow(2*X_DATA[i]*xmp0,2)));
		fac*=units;//change unit to GeV
		//fac = fac /(units * 1.e-3);    //[mikrobarn^-1]

		CS_DATA[i] = fac*CS_DATA[i];
		ERR_DATA[i] = fac*ERR_DATA[i];
		
		if((X_DATA[i]<=X_MAX)&&( Q2_DATA[i]<=Q2_MAX)){
			//fprintf(stdout, "%d: %lE %lE %lE %lE %lE\n",i+1,*(X_DATA+i), *(Y_DATA+i), *(Q2_DATA+i), *(CS_DATA+i), *(ERR_DATA+i));
			i++;
		}	
		j++;
	}
	fprintf(stdout, "%d points, %d\n",i,j);
	fclose(file);
	
	////////////////////////////////////////////////////////////////////////
	
#if NEW_DATA==1 
	N_DATA=i;
	//import_points(X_DATA,Q2_DATA);
	return(0);
#endif	
	printf("load_data:: code not finished...\n");
	return(1);
	//////////////////////////////ZEUS 05//////////////////////////////////
	file =fopen("../data/h1zeus_05.dat","r");
	while((!feof(file) )&&(j<331)){
		fscanf(file, "%lE %lE %lE %lE %lE", (X_DATA+i), (Y_DATA+i) , ( Q2_DATA+i), (CS_DATA+i), (ERR_DATA+i) ) ;
		
		if((X_DATA[i]<=X_MAX)&&( Q2_DATA[i]<=Q2_MAX)){
			//fprintf(stdout, "%d: %lE %lE %lE %lE %lE\n",i+1,*(X_DATA+i), *(Y_DATA+i), *(Q2_DATA+i), *(CS_DATA+i), *(ERR_DATA+i));
			i++;
		}
	j++;	
	}
	fprintf(stdout, "%d points, %d\n",i,j);
	fclose(file);
	

	//////////////////////////H1 LOW 2001////////////////////
	fprintf(stdout, "HERA low 2001\n");
	file = fopen("../data/h1low2001.dat","r");
	j=0;
	while((!feof(file))&&( j<62)){
	
		fscanf(file,"   %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE", &Q2_DATA[i], &X_DATA[i], &Y_DATA[i],&dum,&dum,&CS_DATA[i], &ERR_DATA[i],&dum,&dum,&dum);
            	ERR_DATA[i]= RESCALE*ERR_DATA[i]*CS_DATA[i]/100;
            	CS_DATA[i]= RESCALE*CS_DATA[i];
	
		if((X_DATA[i]<=X_MAX)&&( Q2_DATA[i]<=Q2_MAX)){
			//fprintf(stdout, "%d: %lE %lE %lE %lE %lE\n",i+1,*(X_DATA+i), *(Y_DATA+i), *(Q2_DATA+i), *(CS_DATA+i), *(ERR_DATA+i));
			i++;
		}
			
		j++;
	}
	fprintf(stdout, "%d points, %d\n",i,j);
	fclose(file);
	
	/////////////////////////////H1 2001///////////////////////
	fprintf(stdout, "H1 2001\n");

	file = fopen("../data/h12001.dat","r");  
	j=0;  
	while((!feof(file))&&( j< 71)){
		fscanf(file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &Q2_DATA[i],&X_DATA[i],&Y_DATA[i],&dum,&dum, &CS_DATA[i], &ERR_DATA[i],&dum,&dum,&dum);
		ERR_DATA[i]= RESCALE*ERR_DATA[i]*CS_DATA[i]/100;
		CS_DATA[i]= RESCALE*CS_DATA[i];
		if((X_DATA[i]<=X_MAX)&&( Q2_DATA[i]<=Q2_MAX)){
			//fprintf(stdout, "%d: %lE %lE %lE %lE %lE\n",i+1,*(X_DATA+i), *(Y_DATA+i), *(Q2_DATA+i), *(CS_DATA+i), *(ERR_DATA+i));
			i++;
		}
		j++;	
	}
	fprintf(stdout, "%d points, %d\n",i,j);
	fclose(file);
	
	/////////////////////////////BPT 97///////////////////////
	fprintf(stdout, "BPT 97\n");

	file = fopen("../data/bpt97.dat","r");
	j=0;
	
	while((!feof(file))&&(j<70)){
	    fscanf(file,
		   "%lf %lE %lf %lf %lf %lf %lf", 
		   &Q2_DATA[i],&X_DATA[i],&Y_DATA[i],&CS_DATA[i],
                   &stat,&sysup,&sysdown); 
            dum = maximum(sysup,sysdown);
            ERR_DATA[i] = sqrt(stat*stat+dum*dum);
            if((X_DATA[i]<=X_MAX)&&( Q2_DATA[i]<=Q2_MAX)){
			//fprintf(stdout, "%d: %lE %lE %lE %lE %lE\n",i+1,*(X_DATA+i), *(Y_DATA+i), *(Q2_DATA+i), *(CS_DATA+i), *(ERR_DATA+i));
			i++;
		}
		j++;	
	}
	fprintf(stdout, "%d points, %d\n",i,j);
	fclose(file);
	
	///////////////////////////ZEUS 2001//////////////////////
	fprintf(stdout, "zeus 2001\n");
	file = fopen("../data/zeus2001.dat","r");
	j=0;
	while((!feof(file) )&&(j<242) ){
		fscanf(file,
		   "    %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
		   &Q2_DATA[i],&X_DATA[i],&CS_DATA[i],&dum, 
                   &sysup,&sysdown,&stat,
                   &dum,&dum,&dum,&dum,&dum,&dum,&dum,&dum,&dum,&dum,&dum); 
            stat  = stat*CS_DATA[i]/100;
            sysup = sysup*CS_DATA[i]/100;
            ERR_DATA[i] = sqrt(stat*stat+sysup*sysup);
            //printf("%f  %f \n",Q2_DATA[i], X_DATA[i]);
	if((X_DATA[i]<=X_MAX)&&( Q2_DATA[i]<=Q2_MAX)){
			//fprintf(stdout, "%d: %lE %lE %lE %lE %lE\n",i+1,*(X_DATA+i), *(Y_DATA+i), *(Q2_DATA+i), *(CS_DATA+i), *(ERR_DATA+i));
			i++;
		}	
		j++;
	}
	fprintf(stdout, "%d points, %d\n",i,j);
	fclose(file);
	

	return(0);
}
