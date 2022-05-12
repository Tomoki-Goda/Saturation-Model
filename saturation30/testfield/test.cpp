#include<iostream>
#include"./DIS-cross-section.h"


using namespace std;

int main(){
	cross_section sigma_l;
	sigma_l.show_parameters("extra");

	sigma_l.set_type('l',1);
	sigma_l.show_parameters("other");
	
	return(0);

}
