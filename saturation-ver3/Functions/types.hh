



#if SIGMA_APPROX<0
	//typedef Interpolate_Collinear_Gluon COLGLU;
	typedef Chebyshev_Collinear_Gluon COLGLU;
#else 	
	typedef Collinear_Gluon COLGLU;
#endif
#if MODEL==1
	typedef Sigma<COLGLU> SIGMA ;
#else 
	typedef Sigma SIGMA ;
#endif

#if ((SIGMA_APPROX==0)||(SIGMA_APPROX==-1))//negative means xg is approximated.
//1 grid of sigma
//0 no approx
//-1 only the xg --- default
//-2 BOTH
	#if MODEL==1
		typedef Gluon_Integrand<SIGMA> DSIGMA;
	#else 
		typedef Gluon_Integrand DSIGMA;
	#endif
#elif ((SIGMA_APPROX>0)||(SIGMA_APPROX<-1))
	//typedef Gluon_Integrand<SIGMA> DSIGMA;
	typedef Laplacian_Sigma DSIGMA;
#endif

#if ((GLUON_APPROX==0)&&(MODEL==0))
	typedef Gluon_GBW GLUON;
#else
	typedef Dipole_Gluon<DSIGMA> GLUON;
#endif
