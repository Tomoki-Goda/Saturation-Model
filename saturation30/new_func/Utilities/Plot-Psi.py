#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import math

savedir="/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes"

#dataarray=np.array([]);
def import_array(name ):
	dataarray=[];
	with open(name,"r") as fi:
	#with open("./pwf.txt","r") as fi:
		for line in fi:
			#data=line.readline();
			data=line.strip().split('\t')
			#data=[math.log(float(j))/math.log(10) for j in data]
			#data=[float(j) for j in data ]
			data=[float(data[0]),float(data[1]) ]
			dataarray.append(data)
	dataarray=np.array(dataarray)
	dataarray=np.transpose(dataarray)
	return(dataarray)
	
def import_array_slope(name,xpow ):
	dataarray=[];
	with open(name,"r") as fi:
		for line in fi:
			#data=line.readline();
			data=line.strip().split('\t')
			#data=[math.log(float(j))/math.log(10) for j in data]
			data=[float(j) for j in data]
			data=[data[0],data[2]]
			dataarray.append(data)
	dataarray=np.array(dataarray)
	dataarray=np.transpose(dataarray)
	return(dataarray)
	
#name="NewRun"	
#name="FirstRun"
#name="ThetaOff"
#path= "./"+name+"/"
################################## x dependence ###################################
dataarray=import_array("./photon-psi10.txt")
plt.plot(dataarray[0],dataarray[1],label="Q2=10",ls="-")
dataarray=import_array("./photon-psi50.txt")
plt.plot(dataarray[0],dataarray[1],label="Q2=50",ls="-")
dataarray=import_array("./photon-psi100.txt")
plt.plot(dataarray[0],dataarray[1],label="Q2=100",ls="-")
dataarray=import_array("./photon-psi500.txt")
plt.plot(dataarray[0],dataarray[1],label="Q2=500",ls="-")

plt.xscale('log')
plt.yscale('log')
plt.grid(True)

plt.xlabel("r")
plt.ylabel("|psi|^2")
plt.title("z integrated |Psi|^2")


plt.legend()
plt.savefig(savedir+"/photon-psi.png")
#plt.show()
plt.clf()
		

