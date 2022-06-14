#include<math.h>
#include<stdlib.h>
#include<stdio.h>

#ifndef PI
	#define PI 3.1415
#endif

#ifndef CHEBYSHEV_TEST 
	#define CHEBYSHEV_TEST 0
#endif
////////////////// usage/////////////////
// run
//	cheb_coeff(func, par ,degree, dim ,coeff)
// to produce coeff. 
//then use  
//	chebyshev(degree, dim,  coeff , args )
//////////////////////////////////////////

///////////////////////////Change of variables to -1 1 ///////////////////////
//double change_var_compactify(double min,double max, double val){
double change_var_revert(double min,double max, double val){
	//for val =[-1,1]  return value between min and max
	if((val>1)||(val<-1)|| (min>max)){
		char str[100];
		printf("wrong input for change_var_revert\n val=%f\t [%f, %f] \n",val, min,max);
		printf("%f\t%f\t %f \n type something and enter to continue\n",val, min,max);
		scanf("%s ",str);
		printf("continue");
		//return 0.0;
	}
	return ( (val/2) *(max-min)+ (max+min)/2 );
}
//double change_var_revert(double min,double max, double val){
double change_var_compactify(double min,double max, double val){
	//for val =[min,max]  return value between -1 and 1
	if(min>max){
		char str[100];
		printf("wrong input for change_var_compactify\n val=%f\t [%f, %f] \n",val, min,max);
		printf("%f\t%f\t %f \n type something and enter to continue\n",val, min,max);
		scanf("%s ",str);
		printf("continue");
		return 0.0;
	}
	
	return (2*((val-min)/(max-min)) -1);
	
}
////////////////////////////////log version ////////////////////////////
double change_var_revert_log(double min,double max, double val){
	//for val =[-1,1]  return value between min and max
	double ret=min*pow((max/min),(1.0-val)/2)  ;
	
	return (ret);
}
//double change_var_revert(double min,double max, double val){
double change_var_compactify_log(double min,double max, double val){
	double ret=1-2*(log(val/min) /log(max/min)) ;
	
	return (ret);	
}
//////////////////////////////////////////////////////////////////////////////////
void ind_vec_increment(unsigned * vec, const unsigned * max, unsigned dim){
//increment the vec
// like in the way increasing (hr,min,sec) second by second. 	
	for(int i=(dim-1);i>=0;i--){
		if(*(vec+i) != (*(max+i)-1) ){
		//if it is not max then this will be in creased. if it is then back to zero and next index will be increased.
			(*(vec+i))++;
			break;
		}else{
			(*(vec+i))=0;
		}
	}
}


void chebyshevT(double x,unsigned degree, double * T ){
//iterative  definition of chebyshev polynomial from $T_0(x)$ to $T_{degree-1}(x)$  . 
	*(T)=1;
	*(T+1)=x;
	for(unsigned i=2;i<(degree);i++){
		*(T+i)=2*x*(*(T+i-1))-(*(T+i-2));
	}
}

int kronecker(int i,int j){
	return((i==j)?1:0);
} 


double cheb_c_summand(double func(double * vec,double* par) , double * par, unsigned* ind1, unsigned* ind2, unsigned *degree, unsigned dim ){
//vec, ind1, ind2, degree are vector of length dim.
// sum over ind2 is the coefficients in the chebyshev
	double argvec[dim];
	double factor=1;
	double cosarg; 
	double val;
	
	for(unsigned i=0;i<dim;i++){
		cosarg=PI*(ind2[i]+0.5)/( degree[i] );
		argvec[i]=cos( cosarg ) ;
		factor*=cos(ind1[i] * cosarg );
	}
	
	 val=(*func)(argvec,par)*factor;
	 //printf("sijIJ %f\t %f\n",val,factor);
	 return val;
}


double cheb_c(double func(double * vec,double* par) , double * par, unsigned* ind1,unsigned *degree, unsigned dim ){
// do the sum, produce coefficient c of chebyshev \sum c*T
//vec, ind1,  degree are vector of length dim.
	unsigned ind2[dim];//={0};
		
	double val=0;
	double dif=0;
	unsigned len=1;
	
	for(unsigned i=0;i<dim;i++){
		len*=degree[i];
		*(ind2+i)=0;//initialize indices
	}
	for(unsigned j=0;j<len; j++ ){
//		for(unsigned i=0;i<dim;i++){
//			printf("%d", ind2[i]);
//			(i==(dim-1))?(printf("\n")):(printf(", ")) ;
//		}
		dif=cheb_c_summand( func, par,ind1 ,ind2,degree,dim);
		val+=dif;
#if CHEBYSHEV_TEST ==1
		printf("%f\t",dif);
		
		for(unsigned i=0;i<dim;i++){
			printf("ind1 %d/%d\t",ind1[i],degree[i]);
			printf("ind2 %d/%d\t",ind2[i],degree[i]);
		}
		
		printf("\n");
#endif
		ind_vec_increment(ind2,degree,dim);//increment the ind1; like in the way increasing (hr,min,sec) second by second. 
	}
//	printf("\n");
	return val;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// Main fnctions of this file //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void cheb_coeff(double func(double * vec,double* par) , double * par,unsigned *degree,unsigned dim, double* coeff  ){
	unsigned ind1[dim];//={0};
	unsigned len=1;
	
	for(unsigned i=0;i<dim;i++){
		len*=degree[i];
		*(ind1+i)=0;//initialize indices
	}
	for(unsigned j=0;j<len; j++ ){	
//		for(unsigned i=0;i<dim;i++){
//			printf("%d", ind1[i]);
//			(i==(dim-1))?(printf("\n")):(printf(", ")) ;
//		}
		
		(*(coeff+j))=cheb_c( func, par,ind1 ,degree,dim);
		//printf("coeff\t %f\n", (*(coeff+j)) );
		ind_vec_increment(ind1,degree,dim);//increment the ind1; like in the way increasing (hr,min,sec) second by second. 
		//order in wich coeff is stored is {(0,0,0),(0,0,1),(0,0,2)...,(0,1,0),(0,1,1)...}
		
	}
//	printf("\n");
	
	
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// approximate function ////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double chebyshev(unsigned *degree,unsigned dim, double* coeff ,double* args ){

	//printf("Args\t");
	//for( unsigned i=0;i<dim;i++){
	//	printf("%f\t",args[i]);
	//}
	//printf("\n");
	
	unsigned ind1[dim];
	unsigned len=1;
	unsigned lenT=0;
	double res=0;
	double val;
	
	//double *Tlist[dim];
	
	for(unsigned i=0;i<dim;i++){
		len*=degree[i] ;
		lenT+=degree[i] ;
		*(ind1+i)=0;//initialize indices
	}
	//evaluate chebyshev polynomials Ti(x)
	double Tlist[lenT];
	unsigned posit[10]={0};//max dim=10. This is a list containing where the index move to the next dimension.
	
	for(unsigned i=0;i<dim;i++){
		chebyshevT(args[i], degree[i], Tlist+*(posit+i)); 
		*(posit+i+1)=(*(posit+i)+(*(degree+i)) );
		//printf("%d/%d\n",*(posit+i+1),len);
	}//Now Tlist is a concatenated list of T 
	////////////////////Test T//////////////////////
	FILE* file;
	file=fopen("./testchebT.txt","a");
	for(unsigned i=0;i<dim;i++){
		fprintf(file ,"%f\t", args[i] );
		for(unsigned j=0;j<5;j++){
			fprintf(file ,"%f", Tlist[posit[j]+i ]);
			(j==(5-1))?(fprintf(file,"\n")):(fprintf(file ,"\t")) ;
		}
	}
	fclose(file);
	/////////////////////////////////////////////////////////////////////
	
	for(unsigned j=0;j<len; j++ ){
		//posit=0;//position to start counting in tlist since its joined list.
		val=1;
		for(unsigned i=0;i<dim;i++){
			val *=( ((double)(2-kronecker(ind1[i],0)) )/2);
			val *=Tlist[posit[i]+ind1[i]];
			//printf("%f \t ",Tlist[posit[i]+ind1[i]]);
		}
		//printf("%f\n",val) ;
//		for(unsigned i=0;i<dim;i++){
//			printf("%d", ind1[i]);
//			(i==(dim-1))?(printf("\n")):(printf(", ")) ;
//		}
		//printf("%f\t%f\n", val,(*(coeff+j)) );
		val*=(*(coeff+j));
		//printf("%f \n ",(*(coeff+j)) );
		res+=val;	
		ind_vec_increment(ind1,degree,dim);//increment the ind1; like in the way increasing (hr,min,sec) second by second. 
		//order in wich coeff is stored is {(0,0,0),(0,0,1),(0,0,2)...,(0,1,0),(0,1,1)...}
	}

	for(unsigned i=0;i<dim;i++){
		res*=(2.0/degree[i]);	
	}
	//char ch;
	//scanf("%c",&ch);
	return res;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


