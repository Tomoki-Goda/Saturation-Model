#define PI 3.141592653589793238
#define CA 3.0

#if FLAVOUR==0 
#define NF 3.0
#elif FLAVOUR==1
#define NF 4.0
#else// FLAVOUR==2
#define NF 5.0
#endif


#define LQCD2 0.09 // GeV^2 (300MeV)^2
#define Q0 1.0

#define MASS_L2 0.0196
#define MASS_S2 0.0196
#define MASS_C2 1.96
#define MASS_B2 21.16

//#define NORM 1/(2*pow(PI,3))
#define NORM 0.08282813/2 //i don't know what this is but normalization of psi...
