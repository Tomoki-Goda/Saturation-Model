class Gluon_GBW{
	double sigma_0=0,lambda=0,x_0=0;
	//double Q2=0;
	std::string key;
	
	public:
		explicit Gluon_GBW(){
		}
		void init(const double *par){
			//if(key=="gbw"){
				//printf("parameters set\n");
				sigma_0 =(double)par[0];
				lambda	=(double)par[1];
				x_0	=(double)par[2];
			//}else{
			//	std::cout<<"unknown model: "<<key<<std::endl;
			//}

		}
		~Gluon_GBW(){}
		
	public:
		inline double alpha(const double mu2)const{
			return 4.0/(9.0 *log( ((mu2>2*LQCD2)?(mu2):(2.0*LQCD2))/LQCD2));
		}

	public:
		double operator()(const double x,const double k2,double mu2){
			if(x_0<1.0e-5||x_0>1.0e-3){
				return 0;
			}
			if(lambda<0.05||lambda>0.95){
				return 0;
			}
			double Qs2=pow(x_0/x,lambda);
			double val=3.0/(4*PI*PI)*k2/Qs2*exp(-k2/Qs2);
			if(std::isnan(val)==1){
				return(0);
			}
#if RUN==1
			val*=alpha(mu2);
#endif
			return (sigma_0*val) ;
		}
};
