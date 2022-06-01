#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import math

#dataarray=np.array([]);
def import_array(name ):
	dataarray=[];
	with open(name,"r") as fi:
	#with open("./pwf.txt","r") as fi:
		for line in fi:
			#data=line.readline();
			data=line.strip().split('\t')
			#data=[math.log(float(j))/math.log(10) for j in data]
			data=[float(j) for j in data]
			dataarray.append(data)
	dataarray=np.array(dataarray)
	dataarray=np.transpose(dataarray)
	return(dataarray)


dataarray=import_array("./GBWIntegratedHighAccu/10massiveLCB/plot.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw10massive")

dataarray=import_array("./GBWIntegratedHighAccu/500massiveLCB/plot.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw500massive")

dataarray=import_array("./GBWIntegratedHighAccu/10masslessLCB/plot.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw10massless")

dataarray=import_array("./GBWIntegratedHighAccu/500masslessLCB/plot.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw500massless")

#plt.axis([0.1, 10, 0.001, 1.001])
plt.xscale('log')
plt.yscale('log')
plt.grid(True)

plt.xlabel("ln(r)")
plt.ylabel("ln(sigma/sigma_0)")
plt.title("GBW Dipole Cross-Section")

plt.legend()
plt.show()
		
		
		
		
		
