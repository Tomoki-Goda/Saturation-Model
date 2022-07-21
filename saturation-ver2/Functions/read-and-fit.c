///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                  This is the ugliest file of this project...                                            ///
//        this will read data and compute F2 to compare with data to get chisq.                            ///
//         the most important function here is fcn which should be explained in the MINUIT manual.         ///
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#include"control.h"
#include"control-default.h"
#include"constants.h"




//#include"./read-and-fit.h"
#include"./read-and-fit-cheb.h"

int N_SIMPS=N_SIMPS_R;
int N_CHEB=N_CHEB_R;


///////////////////////////////////////////////////////////////////////////////////////
///////////////  now integrate over r with Simpsons method    /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

//extern int N_CHEB, N_SIMP;


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

double compute_chisq(const double *par){
	
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
//	clock_t time;
//	time=clock();
	double cs[N_DATA];//computed cross-section

	static double psi_arr[(5)*(( N_SIMPS_R>N_CHEB_R)? N_SIMPS_R : N_CHEB_R  )*MAXN];
	
	static int n_cheb,n_simp;
	static double prec;
	//if((n_cheb!=N_CHEB)||((n_simp!=N_SIMP)||prec!=SIGMA_PREC) ){
	//if(((n_simp!=N_SIMPS)||prec!=SIGMA_PREC) ){
	if((n_cheb!=N_CHEB)||( n_simp!=N_SIMPS) || ((prec-SIGMA_PREC)>1.0e-10) ){
		generate_psi_set(psi_arr);
		prec=SIGMA_PREC;
		n_cheb=N_CHEB;
		n_simp=N_SIMPS;
		printf(" Integral : %.2e SIMP %d, CHEB %d \nR_MAX %.2e R_MIN %.2e\n",prec,n_simp,n_cheb,(double)R_MAX,R_MIN);
	}
	
	
	generate_data_set(par, psi_arr, cs);
	
	double chiarr[N_DATA];	
	for(unsigned i=0;i<N_DATA;i++){
		//chisq+=pow( ( cs[i]-CS_DATA[i] )/(ERR_DATA[i]),2);
		//chisq+=pow( ( *(cs+i) - *(CS_DATA+i) )/( *(ERR_DATA+i) ),2);
		chiarr[i]=pow( ( *(cs+i) - *(CS_DATA+i) )/( *(ERR_DATA+i) ),2);
	}
	chisq=KBN_sum(chiarr,N_DATA);
	return(chisq );
}

void fcn(const int *npar, const double grad[], double*fcnval, const double *par,const unsigned *iflag,void (*dum)(void) ){
	
	/////////////////////////////////////////
	//for detail see MINUIT documentatin.
	char outline[200];//entry in the log file 
	clock_t time;
	time=clock();
	
	
	static int counter;
	
	sprintf(outline, "%d %d ",counter++,*iflag);
	log_printf(log_file,outline);
	
	//for(unsigned i=0;i<(*npar);i++){
	for(unsigned i=0;i<N_PAR;i++){
		sprintf(outline, "%.2e, ",*(par+i));
		log_printf(log_file,outline);
	}
	
#if (MODEL==1||MODEL==3)	
	approx_xg(par+1);//generate chebyshev coefficients
#endif
//	if(*iflag==3){
//		double error_array[N_DATA];
//		double cs_array[N_DATA];
//		//double r_step=(R_MAX-R_MIN)/(2*N_SIMPS);
//		simpson_error(SAMPLES,error_array);
//		simpson_sum(SAMPLES,cs_array);	
//		for(int i=0;i<N_DATA;i++){
//			sprintf(outline,"Q2=%.2e x=%.2e, Data=%.3e :>  %.3e\t%.3e\n",Q2_DATA[i],X_DATA[i],CS_DATA[i], cs_array[i],error_array[i]);

//			log_printf(log_file,outline);
//		}
//	}else{	
		*fcnval=compute_chisq(par);
	
		time-=clock();
		
		sprintf(outline,"    %.3e (%.3f), in %.1e sec\n",*fcnval,*fcnval/(N_DATA-N_PAR), -((double)time)/CLOCKS_PER_SEC);
		log_printf(log_file,outline);
//	}
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
	file=fopen("./data/hera_tot.dat","r");
	
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
	
#if (NEW_DATA==1) 
	N_DATA=i;
	
	//import_points(X_DATA,Q2_DATA);
	return(i);
#else
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	printf("load_data:: code not finished...\n");
	return(1);
	//////////////////////////////ZEUS 05//////////////////////////////////
	file =fopen("./data/h1zeus_05.dat","r");
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
	file = fopen("./data/h1low2001.dat","r");
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

	file = fopen("./data/h12001.dat","r");  
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

	file = fopen("./data/bpt97.dat","r");
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
	file = fopen("./data/zeus2001.dat","r");
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
#endif
}
