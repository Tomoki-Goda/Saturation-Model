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
		printf("change_var_revert:: wrong input for change_var_revert\n val=%f\t [%f, %f] \n",val, min,max);
	}
	return ( (val/2) *(max-min)+ (max+min)/2 );
}
//double change_var_revert(double min,double max, double val){
double change_var_compactify(double min,double max, double val){
	//for val =[min,max]  return value between -1 and 1
	if(min>max){
		printf("change_var_compactify:: wrong input for change_var_compactify\n val=%f\t [%f, %f] \n",val, min,max);
	}
	return (2*((val-min)/(max-min)) -1);
	
}
////////////////////////////////log version ////////////////////////////
double change_var_revert_log(double min,double max, double val){
	//for val =[-1,1]  return value between min and max
	if((val>1)||(val<-1)|| (min>max)){
		printf("change_var_revert_log:: wrong input for change_var_revert\n val=%f\t [%f, %f] \n",val, min,max);
	}
	return (min*pow((max/min),(1.0-val)/2) ) ;
}
//double change_var_revert(double min,double max, double val){
double change_var_compactify_log(double min,double max, double val){

	if(min>max){
		printf("change_var_compactify_log:: wrong input for change_var_compactify\n val=%f\t [%f, %f] \n",val, min,max);
	}
	return (1-2*(log(val/min) /log(max/min)) );	
}
////////////////////////////////  indices //////////////////////////////////////////////////
// chebyshev approximation is done in multi dimension but to treat in the same way, tensor are treated as list. 
unsigned convert_index(unsigned * index, unsigned * max, unsigned dim){
	unsigned unit;
	unsigned total=0;
	for(unsigned i=0;i<dim;i++){
		unit=1;
		for(unsigned j=0;j<(dim-i-1);j++){
			unit*= max[dim-j-1];	
		}
		total+=unit*index[i];	
	}
	return(total);
}
/////////////
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

///////////////////////////////////////////////////////////////////////////////////

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

void sample(double func(double *,double* ), double * par,  unsigned * degree, unsigned dim,  double* sample_arr){
///////////////////precompute the function. i.e make the function descrete. /////////////////////////
// degree is a list of how many term you want to go in each variable. dim is dimension( no. of variables). 
//then sample_arr should have length = prod( degree[i] ).
//////////////////////////////////////////////////////////////////////////////////////////////////////
	unsigned index[dim];
	double argvec[dim];
	unsigned max=1;
	double cosarg;
	for(unsigned i=0;i<dim;i++){
		index[i]=0;
		max*=degree[i];
	}
	
	for(unsigned i=0;i<max;i++){
		for(unsigned j=0;j<dim;j++){
			cosarg=PI*(index[j]+0.5)/( degree[j] );
			argvec[j]=cos( cosarg ) ;
		}
		
		(*(sample_arr+i)) =(*func)(argvec,par);
		if(convert_index(index,degree,dim)!=i){
			printf("sample:: Index conversion failed %d?=%d\n", convert_index(index,degree,dim),i);
		}
		ind_vec_increment(index,degree,dim);
	}
}


double cheb_c_summand(double * sample_arr, unsigned* ind1, unsigned* ind2, unsigned *degree, unsigned dim ){
//vec, ind1, ind2, degree are vector of length dim.
// sum over ind2 is the coefficients in the chebyshev
	double factor=1;
	double cosarg; 
	double val;
	
	for(unsigned i=0;i<dim;i++){
		cosarg=PI*(ind2[i]+0.5)/( degree[i] );
		factor*=cos(ind1[i] * cosarg );
	}
	 val=(sample_arr[ convert_index(ind2,degree,dim) ]) *factor;
	 return val;
}


double cheb_c(/*double func(double * vec,double* par) */double * sample_arr,/*double * par,*/ unsigned* ind1,unsigned *degree, unsigned dim ){
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
		dif=cheb_c_summand(sample_arr, ind1 ,ind2,degree,dim);
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
	double sample_arr[len];
	sample(func,  par, degree, dim,   sample_arr);
	
	for(unsigned j=0;j<len; j++ ){	
		(*(coeff+j))=cheb_c(sample_arr,ind1 ,degree,dim);
		ind_vec_increment(ind1,degree,dim);//increment the ind1; like in the way increasing (hr,min,sec) second by second. 
		//order in wich coeff is stored is {(0,0,0),(0,0,1),(0,0,2)...,(0,1,0),(0,1,1)...}
		
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// approximate function ////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double chebyshev(unsigned *degree,unsigned dim, double* coeff ,double* args ){

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
		if(args[i]>1) {
			printf("chebyshev:: arg too large %f at position %d\n",args[i], i);
			(args[i])-=2;			
		}
		if(args[i]<-1) {
			printf("chebyshev:: arg too small %f at position %d\n",args[i], i);
			(args[i])+=2;			
		}
		chebyshevT(args[i], degree[i], Tlist+*(posit+i)); 
		*(posit+i+1)=(*(posit+i)+(*(degree+i)) );
		//printf("%d/%d\n",*(posit+i+1),len);
	}//Now Tlist is a concatenated list of T 
	
	for(unsigned j=0;j<len; j++ ){
		//posit=0;//position to start counting in tlist since its joined list.
		val=1;
		for(unsigned i=0;i<dim;i++){
			val *=( ((double)(2-kronecker(ind1[i],0)) )/2);
			val *=Tlist[posit[i]+ind1[i]];
		}
		val*=(*(coeff+j));
		
		res+=val;	
		ind_vec_increment(ind1,degree,dim);
	}

	for(unsigned i=0;i<dim;i++){
		res*=(2.0/degree[i]);	
	}
	return res;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


