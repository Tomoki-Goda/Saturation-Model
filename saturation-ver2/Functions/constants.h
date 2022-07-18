#define PI 3.141592653589793238
#define CA 3

#define GAMMA_E 0.577215664901532860606512090082402431042


#if FLAVOUR==0 
#define NF 3
#elif FLAVOUR==1
#define NF 4
#else// FLAVOUR==2
#define NF 5
#endif


#define LQCD2 0.09 // GeV^2 (300MeV)^2
#define Q0 1.0


//#define NORM 1/(2*pow(PI,3))
#define NORM 0.041414065//0.08282813//i don't know what this is but normalization of psi...
