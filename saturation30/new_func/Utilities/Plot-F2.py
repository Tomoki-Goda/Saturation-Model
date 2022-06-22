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
	

name="NewRun"	
#name="FirstRun"
#name="ThetaOff"
path= "./"+name+"/"
################################## x dependence ###################################
dataarray=import_array(path+"GBS/Mass0.0-Qup100-Model2-Sud2-rfix0/pointsF2-x4.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 100 massless log10x=-4",ls="-")

dataarray=import_array(path+"GBSPert/Mass0.0-Qup100-Model2-Sud1-rfix0/pointsF2-x4.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 100 massless log10x=-4",ls="-.")

dataarray=import_array(path+"GBWIntegrated/Mass0.0-Qup100-Model2-Sud0-rfixDefault/pointsF2-x4.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 100 massless rfix log10x=-4",ls=":")

dataarray=import_array(path+"BGK/Mass0.0-Qup100-Model1-SudDefault-rfixDefault/pointsF2-x4.txt")
plt.plot(dataarray[0],dataarray[1],label="bgk 100 massless log10x=-4",ls="--")

dataarray=import_array(path+"GBS/Mass0.0-Qup100-Model2-Sud2-rfix0/pointsF2-x6.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 100 massless log10x=-6",ls="-")

dataarray=import_array(path+"GBSPert/Mass0.0-Qup100-Model2-Sud1-rfix0/pointsF2-x6.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 100 massless log10x=-6",ls="-.")

dataarray=import_array(path+"GBWIntegrated/Mass0.0-Qup100-Model2-Sud0-rfixDefault/pointsF2-x6.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 100 massless rfix log10x=-6",ls=":")


dataarray=import_array(path+"BGK/Mass0.0196-Qup100-Model1-SudDefault-rfixDefault/pointsF2-x6.txt")
plt.plot(dataarray[0],dataarray[1],label="bgk pert 100 massive log10x=-6",ls="--")


plt.xscale('log')
plt.yscale('log')
plt.grid(True)

plt.xlabel("Q2")
plt.ylabel("F2")
plt.title("F2 curve")


plt.legend()
plt.savefig("/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots"+name+"/F2100.png")
#plt.show()
plt.clf()
		

################################## slope ###################################
dataarray=import_array_slope(path+"GBS/Mass0.0-Qup100-Model2-Sud2-rfix0/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs 100 massless log10x=-5",ls="-")

dataarray=import_array_slope(path+"GBSPert/Mass0.0-Qup100-Model2-Sud1-rfix0/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs pert 100 massless log10x=-5",ls="-.")

dataarray=import_array_slope(path+"GBWIntegrated/Mass0.0-Qup100-Model2-Sud0-rfixDefault/pointsF2-x5.txt",-5)
plt.plot(dataarray[0],dataarray[1],label="gbw 100 massless rfix log10x=-5",ls=":")

dataarray=import_array_slope(path+"BGK/Mass0.0-Qup100-Model1-SudDefault-rfixDefault/pointsF2-x5.txt",-5)
plt.plot(dataarray[0],dataarray[1],label="bgk 100 massless log10x=-5",ls="--")

dataarray=import_array_slope(path+"GBS/Mass0.0196-Qup100-Model2-Sud2-rfix0/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs 100 massive log10x=-5",ls="-")

dataarray=import_array_slope(path+"GBSPert/Mass0.0196-Qup100-Model2-Sud1-rfix0/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs pert 100 massive log10x=-5",ls="-.")

dataarray=import_array_slope(path+"GBWIntegrated/Mass0.0196-Qup100-Model2-Sud0-rfixDefault/pointsF2-x5.txt",-5)
plt.plot(dataarray[0],dataarray[1],label="gbw 100 massive rfix log10x=-5",ls=":")

dataarray=import_array_slope(path+"BGK/Mass0.0196-Qup100-Model1-SudDefault-rfixDefault/pointsF2-x5.txt",-5)
plt.plot(dataarray[0],dataarray[1],label="bgk 100 massive log10x=-5",ls="--")


plt.xscale('log')
plt.yscale('linear')
plt.grid(True)

plt.xlabel("Q2")
plt.ylabel("slope")
plt.title("Slope")


plt.legend()
plt.savefig("/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots"+name+"/Slope100.png")
#plt.show()
plt.clf()	

################################## slope x dep ###################################
dataarray=import_array_slope(path+"GBS/Mass0.0196-Qup100-Model2-Sud2-rfix0/pointsF2-x4.txt",-4 )
plt.plot(dataarray[0],dataarray[1],label="gbs 100 massive log10x=-4",ls="-")

dataarray=import_array_slope(path+"GBSPert/Mass0.0196-Qup100-Model2-Sud1-rfix0/pointsF2-x4.txt",-4 )
plt.plot(dataarray[0],dataarray[1],label="gbs pert 100 massive log10x=-4",ls="-.")

dataarray=import_array_slope(path+"GBWIntegrated/Mass0.0196-Qup100-Model2-Sud0-rfixDefault/pointsF2-x4.txt",-4)
plt.plot(dataarray[0],dataarray[1],label="gbw 100 massive rfix log10x=-4",ls=":")

dataarray=import_array_slope(path+"BGK/Mass0.0-Qup100-Model1-SudDefault-rfixDefault/pointsF2-x4.txt",-4)
plt.plot(dataarray[0],dataarray[1],label="bgk 100 massless log10x=-4",ls="--")

dataarray=import_array_slope(path+"GBS/Mass0.0196-Qup100-Model2-Sud2-rfix0/pointsF2-x6.txt",-6 )
plt.plot(dataarray[0],dataarray[1],label="gbs 100 massive log10x=-6",ls="-")

dataarray=import_array_slope(path+"GBSPert/Mass0.0196-Qup100-Model2-Sud1-rfix0/pointsF2-x6.txt",-6 )
plt.plot(dataarray[0],dataarray[1],label="gbs pert 100 massive log10x=-6",ls="-.")

dataarray=import_array_slope(path+"GBWIntegrated/Mass0.0196-Qup100-Model2-Sud0-rfixDefault/pointsF2-x6.txt",-6)
plt.plot(dataarray[0],dataarray[1],label="gbw 100 massive rfix log10x=-6",ls=":")

dataarray=import_array_slope(path+"BGK/Mass0.0196-Qup100-Model1-SudDefault-rfixDefault/pointsF2-x6.txt",-6)
plt.plot(dataarray[0],dataarray[1],label="bgk 100 massive log10x=-6",ls="--")


plt.xscale('log')
plt.yscale('linear')
plt.grid(True)

plt.xlabel("Q2")
plt.ylabel("slope")
plt.title("Slope")


plt.legend()
plt.savefig("/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots"+name+"/Slope100xdep.png")
#plt.show()
plt.clf()
################################## x dependence ###################################
dataarray=import_array(path+"GBS/Mass0.0-Qup500-Model2-Sud2-rfix0/pointsF2-x4.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massless log10x=-4",ls="-")

dataarray=import_array(path+"GBSPert/Mass0.0-Qup500-Model2-Sud1-rfix0/pointsF2-x4.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 500 massless log10x=-4",ls="-.")

dataarray=import_array(path+"GBWIntegrated/Mass0.0-Qup500-Model2-Sud0-rfixDefault/pointsF2-x4.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 500 massless rfix log10x=-4",ls=":")

dataarray=import_array(path+"BGK/Mass0.0-Qup500-Model1-SudDefault-rfixDefault/pointsF2-x4.txt")
plt.plot(dataarray[0],dataarray[1],label="bgk 500 massless log10x=-4",ls="--")

dataarray=import_array(path+"GBS/Mass0.0-Qup500-Model2-Sud2-rfix0/pointsF2-x6.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massless log10x=-6",ls="-")

dataarray=import_array(path+"GBSPert/Mass0.0-Qup500-Model2-Sud1-rfix0/pointsF2-x6.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 500 massless log10x=-6",ls="-.")

dataarray=import_array(path+"GBWIntegrated/Mass0.0-Qup500-Model2-Sud0-rfixDefault/pointsF2-x6.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 500 massless rfix log10x=-6",ls=":")

dataarray=import_array(path+"BGK/Mass0.0196-Qup500-Model1-SudDefault-rfixDefault/pointsF2-x6.txt")
plt.plot(dataarray[0],dataarray[1],label="bgk pert 500 massive log10x=-6",ls="--")


plt.xscale('log')
plt.yscale('log')
plt.grid(True)

plt.xlabel("Q2")
plt.ylabel("F2")
plt.title("F2 curve")


plt.legend()
plt.savefig("/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots"+name+"/F2500.png")
#plt.show()
plt.clf()


################################## slope ###################################
dataarray=import_array_slope(path+"GBS/Mass0.0-Qup500-Model2-Sud2-rfix0/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massless log10x=-5",ls="-")

dataarray=import_array_slope(path+"GBSPert/Mass0.0-Qup500-Model2-Sud1-rfix0/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs pert 500 massless log10x=-5",ls="-.")

dataarray=import_array_slope(path+"GBWIntegrated/Mass0.0-Qup500-Model2-Sud0-rfixDefault/pointsF2-x5.txt",-5)
plt.plot(dataarray[0],dataarray[1],label="gbw 500 massless rfix log10x=-5",ls=":")

dataarray=import_array_slope(path+"BGK/Mass0.0-Qup500-Model1-SudDefault-rfixDefault/pointsF2-x5.txt",-5)
plt.plot(dataarray[0],dataarray[1],label="bgk 500 massless log10x=-5",ls="--")

dataarray=import_array_slope(path+"GBS/Mass0.0196-Qup500-Model2-Sud2-rfix0/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massive log10x=-5",ls="-")

dataarray=import_array_slope(path+"GBSPert/Mass0.0196-Qup500-Model2-Sud1-rfix0/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs pert 500 massive log10x=-5",ls="-.")

dataarray=import_array_slope(path+"GBWIntegrated/Mass0.0196-Qup500-Model2-Sud0-rfixDefault/pointsF2-x5.txt",-5)
plt.plot(dataarray[0],dataarray[1],label="gbw 500 massive rfix log10x=-5",ls=":")

dataarray=import_array_slope(path+"BGK/Mass0.0196-Qup500-Model1-SudDefault-rfixDefault/pointsF2-x5.txt",-5)
plt.plot(dataarray[0],dataarray[1],label="bgk 500 massive log10x=-5",ls="--")


plt.xscale('log')
plt.yscale('linear')
plt.grid(True)

plt.xlabel("Q2")
plt.ylabel("slope")
plt.title("Slope")


plt.legend()
plt.savefig("/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots"+name+"/Slope500.png")
#plt.show()
plt.clf()	
################################## slope x dep ###################################
dataarray=import_array_slope(path+"GBS/Mass0.0196-Qup500-Model2-Sud2-rfix0/pointsF2-x4.txt",-4 )
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massive log10x=-4",ls="-")

dataarray=import_array_slope(path+"GBSPert/Mass0.0196-Qup500-Model2-Sud1-rfix0/pointsF2-x4.txt",-4)
plt.plot(dataarray[0],dataarray[1],label="gbs pert 500 massive log10x=-4",ls="-.")

dataarray=import_array_slope(path+"GBWIntegrated/Mass0.0196-Qup500-Model2-Sud0-rfixDefault/pointsF2-x4.txt",-4)
plt.plot(dataarray[0],dataarray[1],label="gbw 500 massive rfix log10x=-4",ls=":")

dataarray=import_array_slope(path+"BGK/Mass0.0-Qup500-Model1-SudDefault-rfixDefault/pointsF2-x4.txt",-4)
plt.plot(dataarray[0],dataarray[1],label="bgk 500 massless log10x=-4",ls="--")

dataarray=import_array_slope(path+"GBS/Mass0.0196-Qup500-Model2-Sud2-rfix0/pointsF2-x6.txt",-6 )
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massive log10x=-6",ls="-")

dataarray=import_array_slope(path+"GBSPert/Mass0.0196-Qup500-Model2-Sud1-rfix0/pointsF2-x6.txt",-6 )
plt.plot(dataarray[0],dataarray[1],label="gbs pert 500 massive log10x=-6",ls="-.")

dataarray=import_array_slope(path+"GBWIntegrated/Mass0.0196-Qup500-Model2-Sud0-rfixDefault/pointsF2-x6.txt",-6)
plt.plot(dataarray[0],dataarray[1],label="gbw 500 massive rfix log10x=-6",ls=":")

dataarray=import_array_slope(path+"BGK/Mass0.0196-Qup500-Model1-SudDefault-rfixDefault/pointsF2-x6.txt",-6)
plt.plot(dataarray[0],dataarray[1],label="bgk 500 massive log10x=-6",ls="--")


plt.xscale('log')
plt.yscale('linear')
plt.grid(True)

plt.xlabel("Q2")
plt.ylabel("slope")
plt.title("Slope")


plt.legend()
plt.savefig("/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots"+name+"/Slope500xdep.png")
#plt.show()
plt.clf()
################################## x dependence ###################################
dataarray=import_array(path+"GBS/Mass0.0-Qup10-Model2-Sud2-rfix0/pointsF2-x4.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 10 massless log10x=-4",ls="-")

dataarray=import_array(path+"GBSPert/Mass0.0-Qup10-Model2-Sud1-rfix0/pointsF2-x4.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 10 massless log10x=-4",ls="-.")

dataarray=import_array(path+"GBWIntegrated/Mass0.0-Qup10-Model2-Sud0-rfixDefault/pointsF2-x4.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 10 massless rfix log10x=-4",ls=":")

dataarray=import_array(path+"BGK/Mass0.0-Qup100-Model1-SudDefault-rfixDefault/pointsF2-x4.txt")
plt.plot(dataarray[0],dataarray[1],label="bgk 100 massless log10x=-4",ls="--")

dataarray=import_array(path+"GBS/Mass0.0-Qup10-Model2-Sud2-rfix0/pointsF2-x6.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 10 massless log10x=-6",ls="-")

dataarray=import_array(path+"GBSPert/Mass0.0-Qup10-Model2-Sud1-rfix0/pointsF2-x6.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 10 massless log10x=-6",ls="-.")

dataarray=import_array(path+"GBWIntegrated/Mass0.0-Qup10-Model2-Sud0-rfixDefault/pointsF2-x6.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 10 massless rfix log10x=-6",ls=":")


dataarray=import_array(path+"BGK/Mass0.0196-Qup10-Model1-SudDefault-rfixDefault/pointsF2-x6.txt")
plt.plot(dataarray[0],dataarray[1],label="bgk pert 10 massive log10x=-6",ls="--")


plt.xscale('log')
plt.yscale('log')
plt.grid(True)

plt.xlabel("Q2")
plt.ylabel("F2")
plt.title("F2 curve")


plt.legend()
plt.savefig("/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots"+name+"/F210.png")
#plt.show()
plt.clf()


################################## slope ###################################
dataarray=import_array_slope(path+"GBS/Mass0.0-Qup10-Model2-Sud2-rfix0/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs 10 massless log10x=-5",ls="-")

dataarray=import_array_slope(path+"GBSPert/Mass0.0-Qup10-Model2-Sud1-rfix0/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs pert 10 massless log10x=-5",ls="-.")

dataarray=import_array_slope(path+"GBWIntegrated/Mass0.0-Qup10-Model2-Sud0-rfixDefault/pointsF2-x5.txt",-5)
plt.plot(dataarray[0],dataarray[1],label="gbw 10 massless rfix log10x=-5",ls=":")

dataarray=import_array_slope(path+"BGK/Mass0.0-Qup10-Model1-SudDefault-rfixDefault/pointsF2-x5.txt",-5)
plt.plot(dataarray[0],dataarray[1],label="bgk 10 massless log10x=-5",ls="--")

dataarray=import_array_slope(path+"GBS/Mass0.0196-Qup10-Model2-Sud2-rfix0/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs 10 massive log10x=-5",ls="-")

dataarray=import_array_slope(path+"GBSPert/Mass0.0196-Qup10-Model2-Sud1-rfix0/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs pert 10 massive log10x=-5",ls="-.")

dataarray=import_array_slope(path+"GBWIntegrated/Mass0.0196-Qup10-Model2-Sud0-rfixDefault/pointsF2-x5.txt",-5)
plt.plot(dataarray[0],dataarray[1],label="gbw 10 massive rfix log10x=-5",ls=":")

dataarray=import_array_slope(path+"BGK/Mass0.0196-Qup10-Model1-SudDefault-rfixDefault/pointsF2-x5.txt",-5)
plt.plot(dataarray[0],dataarray[1],label="bgk 10 massive log10x=-5",ls="--")


plt.xscale('log')
plt.yscale('linear')
plt.grid(True)

plt.xlabel("Q2")
plt.ylabel("slope")
plt.title("Slope")


plt.legend()
plt.savefig("/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots"+name+"/Slope10.png")
#plt.show()
plt.clf()
################################## slope x dep ###################################
dataarray=import_array_slope(path+"GBS/Mass0.0196-Qup10-Model2-Sud2-rfix0/pointsF2-x4.txt",-4 )
plt.plot(dataarray[0],dataarray[1],label="gbs 10 massive log10x=-4",ls="-")

dataarray=import_array_slope(path+"GBSPert/Mass0.0196-Qup10-Model2-Sud1-rfix0/pointsF2-x4.txt",-4 )
plt.plot(dataarray[0],dataarray[1],label="gbs pert 10 massive log10x=-4",ls="-.")

dataarray=import_array_slope(path+"GBWIntegrated/Mass0.0196-Qup10-Model2-Sud0-rfixDefault/pointsF2-x4.txt",-4)
plt.plot(dataarray[0],dataarray[1],label="gbw 10 massive rfix log10x=-4",ls=":")

dataarray=import_array_slope(path+"BGK/Mass0.0-Qup10-Model1-SudDefault-rfixDefault/pointsF2-x4.txt",-4)
plt.plot(dataarray[0],dataarray[1],label="bgk 10 massless log10x=-4",ls="--")

dataarray=import_array_slope(path+"GBS/Mass0.0196-Qup10-Model2-Sud2-rfix0/pointsF2-x6.txt",-6 )
plt.plot(dataarray[0],dataarray[1],label="gbs 10 massive log10x=-6",ls="-")

dataarray=import_array_slope(path+"GBSPert/Mass0.0196-Qup10-Model2-Sud1-rfix0/pointsF2-x6.txt",-6 )
plt.plot(dataarray[0],dataarray[1],label="gbs pert 10 massive log10x=-6",ls="-.")

dataarray=import_array_slope(path+"GBWIntegrated/Mass0.0196-Qup10-Model2-Sud0-rfixDefault/pointsF2-x6.txt",-6)
plt.plot(dataarray[0],dataarray[1],label="gbw 10 massive rfix log10x=-6",ls=":")


dataarray=import_array_slope(path+"BGK/Mass0.0196-Qup10-Model1-SudDefault-rfixDefault/pointsF2-x6.txt",-6)
plt.plot(dataarray[0],dataarray[1],label="bgk 10 massive log10x=-6",ls="--")



plt.xscale('log')
plt.yscale('linear')
plt.grid(True)

plt.xlabel("Q2")
plt.ylabel("slope")
plt.title("Slope")


plt.legend()
plt.savefig("/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots"+name+"/Slope10xdep.png")
#plt.show()
plt.clf()
		
