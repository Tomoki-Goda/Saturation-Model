#ifndef CHEBYSHEV_HH
#define CHEBYSHEV_HH 1
#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include"./Kahn.hh"
#ifndef PI
	#define PI 3.1415926535897932384626433832795
#endif

#ifndef CHEBYSHEV_TEST 
	#define CHEBYSHEV_TEST 0
#endif
////////////////// usage/////////////////
// run
//	cheby PrepareChebyshev(const unsigned *degree,const unsigned dim)
//  degree is a list specifing number of points in each dimension of dim.
//  cheby class is the object
// template<typename T, typename T2>static 
//    void cheb_coeff(cheby & data, T& func, T2 par )
//  to survey function(a,b) a is vector argument, and b vector parameter
//then use  
//	chebyshev(chby, arg) 
//    to get the value
//  
//////////////////////////////////////////

///////////////////////////Change of variables to -1 1 ///////////////////////
//double change_var_compactify(double min,double max, double val){
static double change_var_revert(double min,double max, double val){
	//for val =[-1,1]  return value between min and max
	if((val>1)||(val<-1)|| (min>max)){
		printf("change_var_revert:: wrong input for change_var_revert\n val=%f\t [%f, %f] \n",val, min,max);
		//return (0);
	}
	double res=((val+1)*max - (val-1)*min)/2; 

	return(res);
	//return ( (val/2) *(max-min)+ (max+min)/2 );
}
//double change_var_revert(double min,double max, double val){
static double change_var_compactify(double min,double max, double val){
	//for val =[min,max]  return value between -1 and 1
	if(min>max){
		printf("change_var_compactify:: wrong input for change_var_compactify\n val=%f\t [%f, %f] \n",val, min,max);
		return (0);
	}
	double res=2*((val-min)/(max-min)) -1;
//	printf("%.6e \n",val-change_var_revert(min,max,res)); 
	return(res );
	
}
////////////////////////////////log version ////////////////////////////
static double change_var_revert_log(double min,double max, double val){
	//for val =[-1,1]  return value between min and max
	if((val>1)||(val<-1)|| (min>max)){
		printf("change_var_revert_log:: wrong input for change_var_revert\n val=%f\t [%f, %f] \n",val, min,max);
		//return (0);
	}
	return (min*pow((max/min),(1.0-val)/2) ) ;
}
//double change_var_revert(double min,double max, double val){
static double change_var_compactify_log(double min,double max, double val){

	if(min>max){
		printf("change_var_compactify_log:: wrong input for change_var_compactify\n val=%f\t [%f, %f] \n",val, min,max);
		return (0);
	}
	return (1-2*(log(val/min) /log(max/min)) );	
}
////////////////////////////Frac///////////////////////////////////
//A little special version. 
//use negative value for max to make the max infinity
static double change_var_revert_frac(double min, double max, double val){
	//for val =[-1,1]  return value between min and max

	if((val>1)||(val<-1)){
		printf("change_var_revert_log:: wrong input for change_var_revert\n val=%f \n",val);
	}
	if(max<0){
		return ((1+2*min-val)/(1+val) ) ;
	}else{
		return ((1+2*min-val)*max/(max+2+max*val) ) ;
	}
}
//double change_var_revert(double min,double max, double val){
static double change_var_compactify_frac(double min, double max, double val){
	if(max<0){
		return ((1+2*min-val)/(1+val) );	
	}else{
		return ((max*(1+2*min)-val*(2+max))/(max*(1+val)) );	
	}
}
////////////////////////////////  indices //////////////////////////////////////////////////
// chebyshev approximation is done in multi dimension but to treat in the same way, tensor are treated as list. 
static unsigned convert_index(const unsigned * index, const unsigned * max, const unsigned dim){
	unsigned unit=1;
	unsigned total=0;
	for(unsigned i=0;i<dim;++i){
		total+=unit*index[dim-i-1];	
		unit*=max[dim-i-1];
	}
	return(total);
}
/////////////
static void ind_vec_increment(unsigned * vec, const unsigned * max, const unsigned dim){
//increment the vec
// like in the way increasing (hr,min,sec) second by second. 	
// in the case of time vec should be int vec[3] max={24,60,60}, and dim=3. 
	for(int i =0; i<dim;++i){
		if((++(vec[dim-i-1]))/max[dim-i-1] ==1){
			vec[dim-i-1]=0;
		}else{
			break;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////
/*
void chebyshevT(double x,unsigned degree, double * T ){
//iterative  definition of chebyshev polynomial from $T_0(x)$ to $T_{degree-1}(x)$  . 
	*(T)=1;
	*(T+1)=x;
	for(unsigned i=2;i<(degree);i++){
		*(T+i)=2*x*(*(T+i-1))-(*(T+i-2));
	}
}
*/
static void chebyshevT(const double x, const unsigned degree, double *T){
//iterative  definition of chebyshev polynomial from $T_0(x)$ to $T_{degree-1}(x)$  . 
//with kahn algorithm
	T[0]=1;
	T[1]=x;
	Kahn accum[3]={Kahn_init(2),Kahn_init(2),Kahn_init(2)};
	//double t, c, t0, t1;
	accum[0]+=T[0];
	accum[1]+=T[1];
	
	for(unsigned i=2;i<(degree);++i){
		Kahn_clear(accum[i-3*(i/3)]);
		accum[i-3*(i/3)]+=accum[i-1-3*((i-1)/3)];
		accum[i-3*(i/3)]*=2*x;
		accum[i-3*(i/3)]-=accum[i-2-3*((i-2)/3)];
		T[i]=Kahn_total(accum[i-3*(i/3)]);
		//T[i]=cos(i*acos(x));
	}
	for(int i=0;i<3;++i){
		Kahn_free(accum[i]);
	}
	/*for(unsigned i=2;i<(degree);i++){
		*(T+i)=cos(i*acos(x));
	}*/
}
static void chebyshevU(const double x, const unsigned degree, double * U ){
//iterative  definition of chebyshev polynomial from $T_0(x)$ to $T_{degree-1}(x)$  . 
//with kahn algorithm
	U[0]=1;
	U[1]=2*x;
	Kahn accum[3]={Kahn_init(2),Kahn_init(2),Kahn_init(2)};
	//double t, c, t0, t1;
	accum[0]+=U[0];
	accum[1]+=U[1];
	//double t, c, t0, t1;

	for(unsigned i=2;i<(degree);++i){
		Kahn_clear(accum[i-3*(i/3)]);
		accum[i-3*(i/3)]+=accum[i-1-3*((i-1)/3)];
		accum[i-3*(i/3)]*=2*x;
		accum[i-3*(i/3)]-=accum[i-2-3*((i-2)/3)];
		U[i]=Kahn_total(accum[i-3*(i/3)]);
		/*t1=2*x*(*(U+i-1));
		t0=-(*(U+i-2));

		t=t1+t0;
		if(fabs(t1)>fabs(t0)){
			c=(t1-t)+t0;
		}else{
			c=(t0-t)+t1;
		}
		if(fabs(c)>1.0e-10){
			printf("accum= %.3e %.3e \n",c,t);
		}
		*(U+i)=t+c;*/
		
	}
	for(int i=0;i<3;++i){
		Kahn_free(accum[i]);
	}
}

static void chebyshev_1(const double x, const unsigned degree, double *T1){
	double U[degree];
	chebyshevU(x,degree,U);
	T1[0]=0;
	for(int i=1;i<degree;i++){
		T1[i]=i*U[i-1];
	}
}

static void chebyshev_2(const double x,const unsigned degree, double *T2){
	double U[degree], T[degree];
	chebyshevT(x,degree,T);
	chebyshevU(x,degree,U);
	T2[0]=0;
	T2[1]=0;
	for(int i=2;i<degree;i++){
		T2[i]=i*((i+1)*T[i]-U[i])/(x*x-1);
	}
}






inline int kronecker(const int i,const int j){
	return( ( (i==j)?1:0)  );
} 

template<typename T, typename T2>static void sample(T & func, T2 par,  const unsigned * degree, const unsigned dim,  double* sample_arr){
///////////////////precompute the function. i.e make the function descrete. /////////////////////////
// degree is a list of how many term you want to go in each variable. dim is dimension( no. of variables). 
//then sample_arr should have length = prod( degree[i] ).
//////////////////////////////////////////////////////////////////////////////////////////////////////
	unsigned index[dim];
	//double *argvec;
	unsigned max=1;
	double cosarg;
	for(unsigned i=0;i<dim;i++){
		index[i]=0;
		max*=degree[i];
	}
	//argvec=(double*)malloc(max*dim*sizeof(double));
	double argvec[max*dim];
	for(unsigned i=0;i<max;i++){
		for(unsigned j=0;j<dim;j++){
			cosarg=PI*(index[j]+0.5)/( degree[j] );
			*(argvec+i*dim+j)=cos( cosarg ) ;
		}
		//sample_arr[i]=func(argvec,par);
		if(convert_index(index,degree,dim)!=i){
			printf("sample:: Index conversion failed %d?=%d\n", convert_index(index,degree,dim),i);
		}
		ind_vec_increment(index,degree,dim);
	}

#pragma omp parallel
{
#pragma omp for schedule(dynamic)
	for(unsigned i=0;i<max;i++){
		sample_arr[i]=func(argvec+i*dim,par);
	}
}
	
	//free(argvec);
}


static double cheb_c_summand(const double * sample_arr,const unsigned* ind1, const unsigned* ind2,const unsigned *degree, const unsigned dim ){
//vec, ind1, ind2, degree are vector of length dim.
// sum over ind2 is the coefficients in the chebyshev
	double factor=1;
	double cosarg; 
	double val;
	
	for(unsigned i=0;i<dim;++i){
		cosarg=PI*(ind2[i]+0.5)/( degree[i] );
		factor*=cos(ind1[i] * cosarg );
	}
	 val=(sample_arr[ convert_index(ind2,degree,dim) ]) *factor;
	 return val;
}


static double cheb_c(const double * sample_arr, const unsigned* ind1,const unsigned *degree, const unsigned dim ){
// do the sum, produce coefficient c of chebyshev \sum c*T 
// coefficient c for indices ind1[dim] 
//vec, ind1,  degree are vector of length dim.
	unsigned ind2[dim];//={0};
	double val=0;
	double dif=0;
	unsigned len=1;
	
	for(unsigned i=0;i<dim;i++){
		len*=degree[i];
		*(ind2+i)=0;//initialize indices
	}
	double arr[len];		
	for(unsigned j=0;j<len; j++ ){
		//dif=cheb_c_summand(sample_arr, ind1 ,ind2,degree,dim);
		//val+=dif;
		arr[j]=cheb_c_summand(sample_arr, ind1 ,ind2,degree,dim);

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
	val=Kahn_list_sum(arr,len);
	static int licz=0;
	if(isnan(val)!=0){
		val=0;
		if((licz++)<7){
			printf("cheb_c:: nan encountered. Returning 0.\n");
		}
	}
	return val;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// Main functions of this file //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef struct{
	double* coeff;
	//double *T;
	//double *U;
       	unsigned *degree;
       	//unsigned *index;
       	unsigned dim;
     	int len;
	int lenT;
} cheby;

static cheby PrepareChebyshev(const unsigned *degree,const unsigned dim){
	cheby cont={NULL,NULL,1,1,1};
	cont.degree=(unsigned*)calloc(dim,sizeof(unsigned));
	//cont.index=(unsigned*)calloc(dim,sizeof(unsigned));
	cont.dim=dim;
	cont.len=1;
	cont.lenT=0;
	for(int i=0;i<dim;i++){
		cont.degree[i]=degree[i];
		cont.len*=degree[i];
		cont.lenT+=degree[i];
		//printf("%d \n",coeff_len);
	}
	cont.coeff=(double*)calloc(cont.len,sizeof(double));
	//cont.T=(double*)calloc(cont.lenT,sizeof(double));
	//cont.U=(double*)calloc(cont.lenT,sizeof(double));
	//printf("dim= %d end\n",cont.dim);
	return cont;
}
static inline cheby PrepareChebyshev(const unsigned deg,const unsigned dim){
	unsigned degree[dim];
	for(int i =0;i<dim;++i){
		degree[i]=deg;
	}
	return(PrepareChebyshev(degree,dim));
}
static int FreeChebyshev(cheby & cont){
	free(cont.coeff);
	free(cont.degree);
	//free(cont.index);
	//free(cont.T);
	//free(cont.U);
	return 0;
}

template<typename T, typename T2>static void cheb_coeff(cheby & data, T& func, T2 par ){
	unsigned ind[data.dim];
	for(unsigned i=0;i<data.dim;i++){
		ind[i]=0;//initialize indices
	}
	double sample_arr[data.len];
	sample<T,T2>(func,  par, data.degree, data.dim,   sample_arr);
	for(unsigned j=0;j<data.len; j++ ){	
		data.coeff[j]=cheb_c(sample_arr,ind,data.degree,data.dim);
		ind_vec_increment(ind,data.degree,data.dim);
		//increment the ind1; like in the way increasing (hr,min,sec) second by second. 
		//order in wich coeff is stored is {(0,0,0),(0,0,1),(0,0,2)...,(0,1,0),(0,1,1)...}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// approximate function ////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static int chebyshev_reduce(const cheby & data,cheby & newcheb, const double arg,const int redpos ){
	double Tlist[data.degree[redpos]];
	unsigned ind1[data.dim], ind2[newcheb.dim];
	int j=0;
	for(unsigned i=0;i<data.dim;++i){
		ind1[i]=0;//initialize indices
		if(i!=redpos){
			ind2[j]=0;
			newcheb.degree[j++]=data.degree[i];
		}
	}
	chebyshevT(arg, data.degree[redpos], Tlist); 
	////////////////////////////////
	// Contract c[i,j,...] T[j](x)
	////////////////////////////////
	//double new[new.len]
	Kahn sum=Kahn_init(2);
	for(unsigned j=0;j<newcheb.len; ++j){
		int k=0;
		for(int i=0;i<data.dim;++i){
			(i==redpos)?(ind1[i]=0):(ind1[i]=ind2[k++]);
			//printf("%d\t",ind1[i]);
		}//printf("\n");
		Kahn_clear(sum);
		for(int i=0; i<data.degree[redpos];++i){
			sum+=((i==0)?(0.5):(1))*Tlist[i]*data.coeff[convert_index(ind1,data.degree,data.dim)];
			ind1[redpos]++;		
		}
		newcheb.coeff[convert_index(ind2,newcheb.degree,newcheb.dim)]=Kahn_total(sum)*(2.0/data.degree[redpos]);
		ind_vec_increment(ind2,newcheb.degree,newcheb.dim);		
	}
	Kahn_free(sum);
	return 0;
}


static double chebyshev(const cheby & data, const double* args ){
	double T[data.lenT];
	unsigned ind1[data.dim];
	double res=0;
	double val;
	for(unsigned i=0;i<data.dim;++i){
		ind1[i]=0;//initialize indices
	}
	//unsigned posit[10]={0};
	unsigned posit=0;
/////////////////////Prepare T ////////////////////////////	
	for(unsigned i=0;i<data.dim;++i){
		if(args[i]>1||args[i]<-1) {
			printf("chebyshev:: arg out of range. %f at position %d\n",args[i], i);			
		}
		chebyshevT(args[i], data.degree[i], T+posit); 
		posit+=data.degree[i];
	}//Now T is a concatenated list of Ti(x) 
	//double arr[len];
	Kahn sum=Kahn_init(2);

///////////////////Compute//////////////////////////////
	for(unsigned j=0;j<data.len; ++j ){
		//posit=0;//position to start counting in tlist since its joined list.
		val=1;
	 	posit=0;
		for(unsigned i=0;i<data.dim;++i){
			val *=((ind1[i]==0)?(0.5):(1))*T[posit+ind1[i]];
			posit+=data.degree[i];
		}
		val*=data.coeff[j];
		sum+=val;	
		ind_vec_increment(ind1,data.degree,data.dim);
	}
///////////////////Sum///////////////////////////////
	res=Kahn_total(sum);//Kahn_list_sum(arr,len);
	for(unsigned i=0;i<data.dim;++i){
		res*=(2.0/data.degree[i]);	
	}
	Kahn_free(sum);
	return res;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////// derivative  /////////////////////////////////////////////////////////
///////////////////// del is a list n th derivative for args. [1,0,0] for first derivative wrt first arg of three arguments //////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static double d_chebyshev(const cheby & data, const double* args, int* del  ){
	double T[data.lenT];
	unsigned ind1[data.dim];
	double res=0;
	double val;
	for(unsigned i=0;i<data.dim;++i){
		ind1[i]=0;//initialize indices
	}
	unsigned posit[10]={0};//max dim=10. This is a list containing where the index move to the next dimension.
/////////////////////Prepare T ////////////////////////////	
	for(unsigned i=0;i<data.dim;++i){
		if(args[i]>1||args[i]<-1) {
			printf("chebyshev:: arg out of range. %f at position %d\n",args[i], i);			
		}
		if(del[i]==2){
			chebyshev_2(args[i],data.degree[i],T+posit[i]);
		}else if(del[i]==1){
			chebyshev_1(args[i],data.degree[i],T+posit[i]);
		}else{
			chebyshevT(args[i], data.degree[i],T+posit[i]);
		}
		posit[i+1]=posit[i]+data.degree[i] ;
	}//Now Tlist is a concatenated list of T 
	Kahn sum=Kahn_init(2);
///////////////////Compute//////////////////////////////
	for(unsigned j=0;j<data.len; ++j ){
		//posit=0;//position to start counting in tlist since its joined list.
		val=1;
		for(unsigned i=0;i<data.dim;++i){
			val *=((ind1[i]==0)?(0.5):(1))*T[posit[i]+ind1[i]];
		}
		val*=data.coeff[j];
		sum+=val;	
		ind_vec_increment(ind1,data.degree,data.dim);
	}
///////////////////Sum///////////////////////////////
	res=Kahn_total(sum);
	//res=Kahn_list_sum(arr,len);
	for(unsigned i=0;i<data.dim;++i){
		res*=(2.0/data.degree[i]);	
	}
	Kahn_free(sum);
	return res;
}
inline double d_chebyshev(const cheby & data, const double* args, int ord,int pos ){
	int del[5]={0};
	del[pos]=ord;
	return (d_chebyshev(data, args,del));
}
/*
double d_chebyshev(const cheby & data, const double* args, int* del  ){
	
	const unsigned* degree=data.degree;
	const unsigned& dim=data.dim;
	const double *coeff=data.coeff;
	
	unsigned ind1[dim];
	unsigned len=1;
	unsigned lenT=0;
	double res=0;
	double val;
	
	//double *Tlist[dim];
	
	for(unsigned i=0;i<dim;i++){
		//len*=degree[i] ;
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
			getchar();			
		}
		if(args[i]<-1) {
			printf("chebyshev:: arg too small %f at position %d\n",args[i], i);
			(args[i])+=2;			
			getchar();			
		}
		if(del[i]>2){
			printf("DERIVATIVE HIGHER THAN 2 IS NOT IMPLEMENTED\n");
		}
		if(del[i]==2){
			chebyshev_2(args[i],degree[i],Tlist+*(posit+i));
		}else if(del[i]==1){
			chebyshev_1(args[i],degree[i],Tlist+*(posit+i));
		}else{
			chebyshevT(args[i], degree[i], Tlist+*(posit+i));
		}	
		*(posit+i+1)=(*(posit+i)+(*(degree+i)) );
		//printf("%d/%d\n",*(posit+i+1),len);
	}//Now Tlist is a concatenated list of T 
	double arr[len];
	for(unsigned j=0;j<len; j++ ){
		//posit=0;//position to start counting in tlist since its joined list.
		val=1;
		for(unsigned i=0;i<dim;i++){
			val *=( ((double)(2-kronecker(ind1[i],0)) )/2);
			val *=Tlist[posit[i]+ind1[i]];
		}
		val*=(*(coeff+j));
		
		//res+=val;
		arr[j]=val;	
		ind_vec_increment(ind1,degree,dim);
	}
	res=Kahn_list_sum(arr,len);
	for(unsigned i=0;i<dim;i++){
		res*=(2.0/degree[i]);	
	}
	return res;
}
*/
#endif
