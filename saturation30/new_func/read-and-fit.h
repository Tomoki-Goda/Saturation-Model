//#include<stdio.h>
//#include<math.h>
//#include<time.h>

//#include"./constants.h" 

//#include"./DIS-cross-section.h"


#define MAXN 600
#define RESCALE 1.05

#define NEW_DATA 1 //1 mean only reading

// ////////GLOBAL to this file...////////////////
static double X_DATA[MAXN]={0};
static double Y_DATA[MAXN]={0};
static double wdata[MAXN]={0};
static double Q2_DATA[MAXN]={0};
static double CS_DATA[MAXN]={0};
static double ERR_DATA[MAXN]={0};
static unsigned N_DATA;



extern double sigma_DIS( double ,double, double ,double* );
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
///////////////////////fitting////////////////////////////////////////////////

void dum_func(void){
	//do nothing
}
void fcn(int npar, double grad[], double*fcnval, double *par,unsigned iflag,void (*dum)(void) ){
	/////////////////////////////////////////
	//for detail see MINUIT documentatin.
	//////////////////////////////////////////
	////////////   GLOBALS   /////////////////
	//////////////////////////////////////////
	//double * DATA	experimental data array
	//double * X		bjorken x
	//double * Q2		Q^2
	//double * Y		rapidity
	//unsigned N_DATA	number of data 
	/////////////////////////////////////////
	//printf("fcn \n");
	double chisq=0.0;
	double val;
	clock_t time;
	time=clock();
	for(unsigned i=0;i<N_DATA;i++){
		val=sigma_DIS( *(X_DATA +i) ,*(Q2_DATA+i) ,*(Y_DATA+i) ,par);
#if TEST==2
		printf(" %.3e : %.3e \n",*(CS_DATA+i), val );
#endif	
		chisq+=pow( (val-(*(CS_DATA+i)))/(*(ERR_DATA+i)),2);
	}
	time-=clock();
	
	for(unsigned i=0;i<N_PAR;i++){
		fprintf(stdout, "%.3e, ",*(par+i));
	}
	fprintf(stdout,"%.2e / %d = %.3e, in %.2e sec\n",chisq, N_DATA, chisq/N_DATA, -((double)time)/CLOCKS_PER_SEC);

	*fcnval=chisq;
}

///////////////////////        READ           ///////////////////////////////////////
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
	double units =389.40; //2.56819e-3; //micro-barn to GeV^-2
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
	return(0);
#endif	
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
