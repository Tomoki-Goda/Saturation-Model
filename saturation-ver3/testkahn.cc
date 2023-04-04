#include<iostream>
#include<cmath>
#include"./Functions/Kahn.hh"

int main(){
	int N=7;
	Kahn sum=Kahn_init(2);
	//double arr[N]={1.0e+20,1.03e-35,1,-1.0e-25,-1.01e+20,-1,1.0e+18};
	double arr[N]={1.0e+20,1.03e-35,1,-1.0e-25,-1.0e+20,-1,1.0e-20};
	for(int i=0;i<N;i++){
		sum+=arr[i];		
		Kahn_show(sum);
	}
	Kahn_show(sum);
	std::cout<<"total= "<<Kahn_total(sum)<<std::endl;
	for(int i=0;i<N;i++){
		sum+=arr[N-1-i];		
		Kahn_show(sum);
	}
	Kahn_show(sum);
	std::cout<<"total= "<<Kahn_total(sum)<<std::endl;
	for(int i=0;i<N;i++){
		sum+=-arr[i];		
		Kahn_show(sum);
	}
	Kahn_show(sum);
	std::cout<<"total= "<<Kahn_total(sum)<<std::endl;

	return 0;
}
