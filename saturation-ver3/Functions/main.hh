#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<vector>
#include<ctime>
#include<chrono>
#include"./control.h"
#include"./control-default.h"
#include"./constants.h"
#include"./Parameters.hh"

double INT_PREC=1;
int N_APPROX=N_CHEB_R;
#include"kt-formula.hh"
#include"fcn.h"

int Fix_Release(ROOT::Minuit2::FunctionMinimum & min, KtFCN theFCN, int N, int str ,int goal);
