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
	
	
#name="FirstRun"
name="ThetaOff"
path= "./"+name+"/"
################################## x dependence ###################################
dataarray=import_array(path+"GBS/100masslessLCB/pointsF2-x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 100 massless log10x=-3",ls="-")

dataarray=import_array(path+"GBSPert/100masslessLCB/pointsF2-x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 100 massless log10x=-3",ls="-.")

dataarray=import_array(path+"GBWIntegrated/100masslessLCB/pointsF2-x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 100 massless rfix log10x=-3",ls=":")

dataarray=import_array(path+"GBS/100masslessLCB/pointsF2-x5.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 100 massless log10x=-5",ls="-")

dataarray=import_array(path+"GBSPert/100masslessLCB/pointsF2-x5.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 100 massless log10x=-5",ls="-.")

dataarray=import_array(path+"GBWIntegrated/100masslessLCB/pointsF2-x5.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 100 massless rfix log10x=-5",ls=":")


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
dataarray=import_array_slope(path+"GBS/100masslessLCB/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs 100 massless log10x=-5",ls="-")

dataarray=import_array_slope(path+"GBSPert/100masslessLCB/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs pert 100 massless log10x=-5",ls="-.")

dataarray=import_array_slope(path+"GBWIntegrated/100masslessLCB/pointsF2-x5.txt",-5)
plt.plot(dataarray[0],dataarray[1],label="gbw 100 massless rfix log10x=-5",ls=":")


dataarray=import_array_slope(path+"GBS/100massiveLCB/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs 100 massive log10x=-5",ls="-")

dataarray=import_array_slope(path+"GBSPert/100massiveLCB/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs pert 100 massive log10x=-5",ls="-.")

dataarray=import_array_slope(path+"GBWIntegrated/100massiveLCB/pointsF2-x5.txt",-5)
plt.plot(dataarray[0],dataarray[1],label="gbw 100 massive rfix log10x=-5",ls=":")


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
################################## x dependence ###################################
dataarray=import_array(path+"GBS/500masslessLCB/pointsF2-x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massless log10x=-3",ls="-")

dataarray=import_array(path+"GBSPert/500masslessLCB/pointsF2-x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 500 massless log10x=-3",ls="-.")

dataarray=import_array(path+"GBWIntegrated/500masslessLCB/pointsF2-x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 500 massless rfix log10x=-3",ls=":")

dataarray=import_array(path+"GBS/500masslessLCB/pointsF2-x5.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massless log10x=-5",ls="-")

dataarray=import_array(path+"GBSPert/500masslessLCB/pointsF2-x5.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 500 massless log10x=-5",ls="-.")

dataarray=import_array(path+"GBWIntegrated/500masslessLCB/pointsF2-x5.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 500 massless rfix log10x=-5",ls=":")


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
dataarray=import_array_slope(path+"GBS/500masslessLCB/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massless log10x=-5",ls="-")

dataarray=import_array_slope(path+"GBSPert/500masslessLCB/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs pert 500 massless log10x=-5",ls="-.")

dataarray=import_array_slope(path+"GBWIntegrated/500masslessLCB/pointsF2-x5.txt",-5)
plt.plot(dataarray[0],dataarray[1],label="gbw 500 massless rfix log10x=-5",ls=":")

dataarray=import_array_slope(path+"GBS/500massiveLCB/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massive log10x=-5",ls="-")

dataarray=import_array_slope(path+"GBSPert/500massiveLCB/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs pert 500 massive log10x=-5",ls="-.")

dataarray=import_array_slope(path+"GBWIntegrated/500massiveLCB/pointsF2-x5.txt",-5)
plt.plot(dataarray[0],dataarray[1],label="gbw 500 massive rfix log10x=-5",ls=":")



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

################################## x dependence ###################################
dataarray=import_array(path+"GBS/10masslessLCB/pointsF2-x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 10 massless log10x=-3",ls="-")

dataarray=import_array(path+"GBSPert/10masslessLCB/pointsF2-x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 10 massless log10x=-3",ls="-.")

dataarray=import_array(path+"GBWIntegrated/10masslessLCB/pointsF2-x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 10 massless rfix log10x=-3",ls=":")

dataarray=import_array(path+"GBS/10masslessLCB/pointsF2-x5.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 10 massless log10x=-5",ls="-")

dataarray=import_array(path+"GBSPert/10masslessLCB/pointsF2-x5.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 10 massless log10x=-5",ls="-.")

dataarray=import_array(path+"GBWIntegrated/10masslessLCB/pointsF2-x5.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 10 massless rfix log10x=-5",ls=":")


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
dataarray=import_array_slope(path+"GBS/10masslessLCB/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs 10 massless log10x=-5",ls="-")

dataarray=import_array_slope(path+"GBSPert/10masslessLCB/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs pert 10 massless log10x=-5",ls="-.")

dataarray=import_array_slope(path+"GBWIntegrated/10masslessLCB/pointsF2-x5.txt",-5)
plt.plot(dataarray[0],dataarray[1],label="gbw 10 massless rfix log10x=-5",ls=":")

dataarray=import_array_slope(path+"GBS/10massiveLCB/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs 10 massive log10x=-5",ls="-")

dataarray=import_array_slope(path+"GBSPert/10massiveLCB/pointsF2-x5.txt",-5 )
plt.plot(dataarray[0],dataarray[1],label="gbs pert 10 massive log10x=-5",ls="-.")

dataarray=import_array_slope(path+"GBWIntegrated/10massiveLCB/pointsF2-x5.txt",-5)
plt.plot(dataarray[0],dataarray[1],label="gbw 10 massive rfix log10x=-5",ls=":")


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

		
