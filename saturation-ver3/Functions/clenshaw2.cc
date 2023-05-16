#include"clenshaw2.hh"
int  Clenshaw_Curtis::cc_init(const int N){
	this->N=N;
	Kahn vec[N/4+1];
	double t[N/4+1];
	double c[N/2+1];
	if((N/8)*8!=N || N>256){
		printf("N=%d has to be multiple of 8, <= 128\n",N );
	}
	x[0]=1;
	x[1]=cos(PI/N);
	for(int i =2;i<N/2;i++){
		x[i]=cos(i*PI/N);
	}
	x[N/2]=0;

	//double *__restrict__ c=(double*)calloc((N/2+1),sizeof(double));
	c[0]=1;
	c[N/2]=1.0/(1-N*N);
	for(int i=1;i<N/2;i++){
		c[i]=2.0/(1-4*i*i);
	}
	//double *__restrict__ t=(double*)calloc((N/4+1),sizeof(double));
	//double accum[3]={0};
	Kahn accum=Kahn_init(4);

	for(int j=0;j<N/4+1;j++){
		if(j==0){
			t[j]=0;
		}else if(j==N/4){
			t[j]=-1;
		}else{
			t[j]=-sin(2*j*PI/N);
		}
		//printf("%.3e\t",t[j]);
		//t[i*(N/2+1)+j]=((i*j==0)?(1):(cos(2*i*j*PI/N)) );
	}//printf("\n");
	int pos;
	//Kahn * vec=(Kahn*)malloc((N/4+1)*sizeof(Kahn));
	for(int j=0;j<N/4+1;j++){
		vec[j]=Kahn_init(4);
		Kahn_clear(vec[j]);
	}

	for(int i=0;i<N/2+1;i++){
		Kahn_clear(accum);
		wfull[i]=0;
		for(int j=0;j<N/2+1;j++){
			pos=t_pos(i,j,N);
			//accum+= sign(pos) * t[abs(pos)] *c[j];
			vec[abs(pos)]+=sign(pos)*c[j];
		//	printf("%.2e\t",t[abs(pos)]);
		}//printf("\n");
		for(int j=0;j<N/4+1;j++){
			vec[j]*=t[j];
			accum+=vec[j];
			Kahn_clear(vec[j]);
		}
		
		wfull[i]=Kahn_total(accum);
		if(i==0){
			wfull[i]/=2;
		}
		//printf("%.10e\t",wfull[i]);
	}//printf("\n");

	//t[0]=1;
	for(int i=0;i<N/4+1;i++){
		Kahn_clear(accum);
		whalf[i]=0;
		for(int j=0;j<N/4+1;j++){
			//t[j]=cos(i*j*4*PI/N);
			pos=t_pos(i,j,N/2);
			vec[2*abs(pos)]+=sign(pos)*c[j];
			//accum+=sign(pos)*t[2*abs(pos)]*c[j];
		}
		for(int j=0;j<N/8+1;j++){
			vec[2*j]*=t[2*j];
			accum+=vec[2*j];
			Kahn_clear(vec[2*j]);
		}
		whalf[i]=Kahn_total(accum);
		
		//whalf[i]+=accum;
		if(i==0){
			whalf[i]/=2;
		}
		//printf("%.10e\t",whalf[i]);
	}//printf("\n");

	//free(c);
	//free(t);
	Kahn_free(accum);
	for(int j=0;j<N/4+1;j++){
		Kahn_free(vec[j]);
	}
	//free(vec);
	//getchar();

	Kahn sum1=Kahn_init(2),sum2=Kahn_init(2);

	for(int i=0;i<N/2;i++){
		sum1+=2*wfull[i];
		if((i/2)*2==i){
			sum2+=2*whalf[i/2];
		}
	}
	sum1+=wfull[N/2];
	sum2+=whalf[N/4];
	double sum11=Kahn_total(sum1),sum22=Kahn_total(sum2);

	sum2*=-2;
	sum2+=sum1;

	if( fabs(Kahn_total(sum2)/(std::min(fabs(sum11),fabs(sum22))))>1.0e-15){
	//if( (sum1-2*sum2)!=0.0){
		printf("Factor array may be problematic, N = %d, %.3e - %.3e = %.3e\n",N,sum11,2*sum22,Kahn_total(sum2) );
		for(int i=0;i<N/2;i++){
			printf("%.3e ",wfull[i]);
			if((i/2)*2==i){
				printf("%.3e ",2*whalf[i/2]);
			}
			printf("\n");
		}
		getchar();
	}	
	Kahn_free(sum1);
	Kahn_free(sum2);
	//getchar();
	return 0;
}

int Clenshaw_Curtis::cc_init(const int N, const std::string& name,int max_rec){
			cc_init(N);
			this->tag=name;
			this->max_rec=max_rec;
			return 0;
		};
int Clenshaw_Curtis::cc_init(const int N, const std::string& name,int InitDiv,int max_rec){
			cc_init( N,name,max_rec);
			this->InitDiv=InitDiv;
			return 0;
		};

double Clenshaw_Curtis::fixed_cc(const void * par, const double smin,const double smax, Kahn &full, Kahn &half)const{
	const double (&x16)[129]=x;
	const double (&w16)[129]=wfull;
	const double (&w8)[65]=whalf;
	double f[257];
	double arg;
	double scale=(smax-smin)/2;
	double mid=(smax+smin)/2;
	for(int i=0;i<N/2;i++){
		arg=scale*x16[i];
		f[2*i]=cc_integrand(mid-arg,par);
		f[2*i+1]=cc_integrand(mid+arg,par);
	}

	f[N]=cc_integrand(mid,par);
	Kahn_clear(full);
	Kahn_clear(half);
	for(int i=0;i<N/2;i++){
		full+=f[2*i]*w16[i];
		full+=f[2*i+1]*w16[i];
	}
	full+=f[N]*w16[N/2];
	for(int i=0;i<N/4;i++){
		half+=f[4*i]*w8[i];
		half+=f[4*i+1]*w8[i];
	}
	half+=f[N]*w8[N/4];
	full*=2*scale/N;
	half*=4*scale/N;
	return Kahn_total(full);
}

double Clenshaw_Curtis::cc_integrate(const void * par,const double a, const double b, const double eps, const double Aeps)const{
	if(a-b==0.0){
		return(0);
	}
	int MAX_RECURSION=max_rec;
	//double integrate(double(&integrand)(const double*, const void*),const void* args,const double a,const double b,const double eps){
	double sign, max,min;

	if(b>a){
		min=a;
		max=b;
		sign=1;
	}else if(a>b){
		max=a;
		min=b;
		sign=-1;
	}else{
		return(0);
	}
	double smin,smax;
	double valfull,valhalf;
	double scale;
	double arg;
	double total=0;//,total2=0;
	double arr[257];
	Kahn accum=Kahn_init(3);
	Kahn accumfull=Kahn_init(3);
	Kahn accumhalf=Kahn_init(3);//, accum2=Kahn_init(3);
	double increase;
	smin=min;
	smax=(min-min/InitDiv)+max/InitDiv;//init_div;
	int licz=0,licztot=0 , counter=0;
	double error_ratio;//,derivs0=1,derivs1=1;	
	while(1){
	if(counter==MAX_RECURSION){
		smax=smin+2*scale;
		printf("Clenshaw_Curtis:: in \"%s\", evaluated %d times.\n",(tag).c_str(),counter );
		printf("sector size = %.3e\n [%.3e, %.3e] of [%.3e, %.3e] after %d / %d /%d\n",smax-smin,smin,smax,min,max, licz,licztot,MAX_RECURSION);
		//getchar();
		//if(counter>=MAX_RECURSION+2){
			goto Error;
		//}
	}
	++counter;
	++licztot;
	scale=(smax-smin)/2;
	
	fixed_cc(par,smin,smax,accumfull,accumhalf);
	valfull=Kahn_total(accumfull);
	valhalf=Kahn_total(accumhalf);
#if integrate_HH==1		
	if(not(std::isfinite(valfull)&&std::isfinite(valhalf))){
		printf("Clenshaw_Curtis:: in \"%s\" %.3e  %.3e encountered\n",(tag).c_str(),valfull,valhalf);
		goto Error;
	}
#endif
	//accumhalf*=-1;
	//accumhalf+=accumfull;
	error_ratio=fabs( (valfull-valhalf)/(eps*valfull) );
	//if(( fabs(valfull-valhalf)<eps*(fabs(valfull)) ) || (  fabs(valfull-valhalf)< fabs(smax-smin)*Aeps ) ){
	if((error_ratio<1 ) || (  fabs(valhalf-valfull)< Aeps ) ){
		accum+=accumfull;
		++licz;
		counter=0;
		if(fabs(smax-max)==0.0){
			total=sign*Kahn_total(accum);
			Kahn_free(accum);
			Kahn_free(accumfull);
			Kahn_free(accumhalf);
			return(total );
		}
		smin=smax;
		increase=(4*scale);
		smax=((max-(smin+increase)<(increase/2))?(max):(smin+increase));
	}else{
		//IMPROVE
		smax=smin+(scale/(1+pow(2*counter,2)*pow(error_ratio, 0.25 )));
		//smax=smin+(scale/2);
	}

	if(((max-min)-(smax-smin))==(max-min)){
		printf("Clenshaw_Curtis:: in \"%s\", division exceeds limitation. in the domain [%.3e, %.3e] of [%.3e, %.3e] after %d / %d \n",(tag).c_str(), smin,smax,min,max, licz,licztot);						
		printf("valfull= %.3e , valhalf= %.3e \n",valfull,valhalf);
		//getchar();
		goto Error;
	}
	
	}

	Error:
		printf("valfull= %.3e , valhalf= %.3e  diff=%.3e\n",valfull,valhalf,valfull-valhalf);
		const double *x16=x;
		const double *w16=wfull;
		const double *w8=whalf;
		scale=(smax-smin)/2;
		double mid=(smax+smin)/2;
		for(int i=0;i<N/2;i++){
			arg=mid-scale*x16[i];
			printf("%.3e\t%.3e\n",arg,cc_integrand(arg,par));
		}
		arg=mid;
		printf("%.3e\t%.3e\n",arg,cc_integrand(arg,par));
		for(int i=0;i<N/2;i++){
			arg=mid+scale*x16[N/2-i-1];
			printf("%.3e\t%.3e\n",arg,cc_integrand(arg,par));
		}
		printf("\n");
		Kahn_free(accum);
		Kahn_free(accumfull);
		Kahn_free(accumhalf);
		exit(1);
		//getchar();
		return 0;

}

