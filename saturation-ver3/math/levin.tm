#include<wstp.h>
#include"Levin.hh"

Levin levin(50);

inline double LevinAdd(double x){
	int pos=levin.add_term(x);
	return (levin.sum(pos));
}
inline double LevinAccel(int i, int j){
	return (levin.accel(i,j));
}
inline int LevinReset(){
	return (levin.reset());
}
inline double LevinSum(int i){
	return (levin.sum(i));
}
int main(int argc,char** argv){
	return WSMain(argc, argv);
}


:Begin:
:Function: LevinAdd
:Pattern: LevinAdd[x_?NumberQ]
:Arguments: {N[x]}
:ArgumentTypes: {Real}
:ReturnType: Real
:End:

:Begin:
:Function: LevinAccel
:Pattern: LevinAccel[i_?NumberQ,j_?NumberQ]
:Arguments: {N[i],N[j]}
:ArgumentTypes: {Integer, Integer}
:ReturnType: Real
:End:

:Begin:
:Function: LevinReset
:Pattern: LevinReset[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Integer
:End:

:Begin:
:Function: LevinSum
:Pattern: LevinSum[j_?NumberQ]
:Arguments: {N[j]}
:ArgumentTypes: {Integer}
:ReturnType: Real
:End:
