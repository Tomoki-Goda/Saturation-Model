#include<math.h>
#include<stdlib.h>
#include<stdio.h>

#ifndef PI
	#define PI 3.1415
#endif

////////////////// usage/////////////////
// run
//	cheb_coeff(func, par ,degree, dim ,coeff)
// to produce coeff. 
//then use  
//	chebyshev(degree, dim,  coeff , args )
//////////////////////////////////////////


void ind_vec_increment(unsigned * vec, const unsigned * max, unsigned dim){
//increment the vec
// like in the way increasing (hr,min,sec) second by second. 
	for(unsigned i=0;i<dim;i++){
		if((*(vec+i)) != (*(max+i))+1 ){
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


double cheb_c_summand(double func(double * vec,double* par) , double * par, unsigned* ind1, unsigned* ind2, const unsigned *degree, unsigned dim ){
//vec, ind1, ind2, degree are vector of length dim.
// sum over ind2 is the coefficients in the chebyshev
	double argvec[dim];
	double factor=1;
	double cosarg; 
	
	for(unsigned i=0;i<dim;i++){
	cosarg=PI*(ind2[i]+0.5)/( degree[i] );
		argvec[i]=cos( cosarg ) ;
		factor*=cos(ind1[i] * cosarg );
	}
	return( (*func)(argvec,par)*factor);
}


double cheb_c(double func(double * vec,double* par) , double * par, unsigned* ind1,const unsigned *degree, unsigned dim ){
// do the sum, produce coefficient c of chebyshev \sum c*T
//vec, ind1,  degree are vector of length dim.
	unsigned ind2[dim];//={0};
		
	double val=0;
	unsigned len=1;
	
	for(unsigned i=0;i<dim;i++){
		len*=degree[i];
		*(ind2+i)=0;//initialize indices
	}
	for(unsigned j=0;j<len; j++ ){
		val+=cheb_c_summand( func, par,ind1 ,ind2,degree,dim);
		ind_vec_increment(ind2,degree,dim);//increment the ind1; like in the way increasing (hr,min,sec) second by second. 
	}
	return val;
}

/////////////////////////////////////////////// Main fnctions of this file //////////////////////////////////////////////////////////////
void cheb_coeff(double func(double * vec,double* par) , double * par,const unsigned *degree,unsigned dim, double* coeff  ){
	unsigned ind1[dim];//={0};
	unsigned len=1;
	
	for(unsigned i=0;i<dim;i++){
		len*=degree[i];
		*(ind1+i)=0;//initialize indices
	}
	for(unsigned j=0;j<len; j++ ){
		(*(coeff+j))=cheb_c( func, par,ind1 ,degree,dim);
		ind_vec_increment(ind1,degree,dim);//increment the ind1; like in the way increasing (hr,min,sec) second by second. 
		//order in wich coeff is stored is {(0,0,0),(0,0,1),(0,0,2)...,(0,1,0),(0,1,1)...}
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
double chebyshev(const unsigned *degree,unsigned dim, double* coeff ,double* args ){
	unsigned ind1[dim];//0};
	unsigned len=1;
	double res=0;
	double val;
	
	for(unsigned i=0;i<dim;i++){
		len*=degree[i];
		*(ind1+i)=0;//initialize indices
	}
	//evaluate chebyshev polynomials Ti(x)
	double Tlist[len];
	unsigned posit=0;
	
	for(unsigned i=0;i<dim;i++){
		chebyshevT(args[i], degree[i], Tlist+posit);
		posit+=degree[i];
	}//Now Tlist is a concatenated list of T 
	
	for(unsigned j=0;j<len; j++ ){
		posit=0;//position to start counting in tlist since its joined list.
		val=1;
		for(unsigned i=0;i<dim;i++){
			val *=((double)(2-kronecker(ind1[i],0))/2);
			val *=Tlist[posit+ind1[i]];
			posit+=degree[i];
		}
		val*=(*(coeff+j));
		res+=val;	
		ind_vec_increment(ind1,degree,dim);//increment the ind1; like in the way increasing (hr,min,sec) second by second. 
		//order in wich coeff is stored is {(0,0,0),(0,0,1),(0,0,2)...,(0,1,0),(0,1,1)...}
	}
	return res;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


