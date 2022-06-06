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
	
	
name="ThetaOff"
path= "./"+name+"/"
################################## x dependence ###################################
dataarray=import_array(path+"GBSrfix/500masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massless rfix log10x=-3",ls="-")

dataarray=import_array(path+"GBSrfix/500masslessLCB/pointsQ2500x5.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massless rfix log10x=-5",ls="-.")

dataarray=import_array(path+"GBSrfix/10masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 10 massless rfix log10x=-3",ls=":")

dataarray=import_array(path+"GBSrfix/10masslessLCB/pointsQ2500x5.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 10 massless rfix log10x=-5",ls="--")


plt.xscale('log')
plt.yscale('log')
plt.grid(True)

plt.xlabel("r /GeV")
plt.ylabel("sigma/sigma_0")
plt.title("GBS (r=0.5,m_l=0) Dipole Cross-Section Q^2=100GeV")


plt.legend()
plt.savefig("/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots"+name+"/GBSrfix x dependence.png")
#plt.show()
plt.clf()
############################## Q dependence ####################################
dataarray=import_array(path+"GBSrfix/500masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massless rfix Q2=500",ls="-")

dataarray=import_array(path+"GBSrfix/500masslessLCB/pointsQ250x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massless rfix Q2=50",ls="-.")

dataarray=import_array(path+"GBSrfix/50masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 50 massless rfix Q2=500",ls=":")

dataarray=import_array(path+"GBSrfix/50masslessLCB/pointsQ250x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 50 massless rfix Q2=50",ls="--")


plt.xscale('log')
plt.yscale('log')
plt.grid(True)

plt.xlabel("r /GeV")
plt.ylabel("sigma/sigma_0")
plt.title("GBS (r=0.5,m_l=0) Dipole Cross-Section x=0.001")


plt.legend()
plt.savefig("/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots"+name+"/GBSrfix Q2 dependence.png")
#plt.show()
plt.clf()	

############################## Model dependence ####################################
dataarray=import_array(path+"GBS/500masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massless",ls="-")


dataarray=import_array(path+"GBSrfix/500masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massless rmax=0.5",ls="--")

dataarray=import_array(path+"GBSPert/500masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 500 massless",ls="--")

dataarray=import_array(path+"GBSPertrfix/500masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 500 massless rmax=0.5",ls="-.")

dataarray=import_array(path+"GBWIntegrated/500masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 500 massless",ls=":")

plt.xscale('log')
plt.yscale('log')
plt.grid(True)

plt.xlabel("r /GeV")
plt.ylabel("sigma/sigma_0")
plt.title("Dipole Cross-Section Q2=500,  x=0.001")


plt.legend()
plt.savefig("/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots"+name+"/Model dependence 500 massless.png")
#plt.show()
plt.clf()

#########################################massive 500 #######################################################
dataarray=import_array(path+"GBS/500massiveLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massive",ls="-")

dataarray=import_array(path+"GBSrfix/500massiveLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massive rmax=0.5",ls="--")

dataarray=import_array(path+"GBSPert/500massiveLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 500 massive",ls="--")

dataarray=import_array(path+"GBSPertrfix/500massiveLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 500 massive rmax=0.5",ls="-.")

dataarray=import_array(path+"GBWIntegrated/500massiveLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 500 massive",ls=":")


plt.xscale('log')
plt.yscale('log')
plt.grid(True)

plt.xlabel("r /GeV")
plt.ylabel("sigma/sigma_0")
plt.title("Dipole Cross-Section Q2=500,  x=0.001")


plt.legend()
plt.savefig("/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots"+name+"/Model dependence 500 massive.png")
#plt.show()
plt.clf()

############################################ massless 50 #################################################################	
dataarray=import_array(path+"GBS/50masslessLCB/pointsQ250x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 50 massless",ls="-")

dataarray=import_array(path+"GBSrfix/50masslessLCB/pointsQ250x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 50 massless rmax=0.5",ls="--")

dataarray=import_array(path+"GBSPert/50masslessLCB/pointsQ250x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 50 massless",ls="--")

dataarray=import_array(path+"GBSPertrfix/50masslessLCB/pointsQ250x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 50 massless rmax=0.5",ls="-.")

dataarray=import_array(path+"GBWIntegrated/50masslessLCB/pointsQ250x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 50 massless",ls=":")

plt.xscale('log')
plt.yscale('log')
plt.grid(True)

plt.xlabel("r /GeV")
plt.ylabel("sigma/sigma_0")
plt.title("Dipole Cross-Section Q2=50,  x=0.001")


plt.legend()
plt.savefig("/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots"+name+"/Model dependence 50 massless.png")
#plt.show()
plt.clf()

############################################ massive 50 ####################################################################
dataarray=import_array(path+"GBS/50massiveLCB/pointsQ250x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 50 massive",ls="-")

dataarray=import_array(path+"GBSrfix/50massiveLCB/pointsQ250x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 50 massive rmax=0.5",ls="--")

dataarray=import_array(path+"GBSPert/50massiveLCB/pointsQ250x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 50 massive",ls="--")

dataarray=import_array(path+"GBSPertrfix/50massiveLCB/pointsQ250x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 50 massive rmax=0.5",ls="-.")

dataarray=import_array(path+"GBWIntegrated/50massiveLCB/pointsQ250x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 50 massive",ls=":")


plt.xscale('log')
plt.yscale('log')
plt.grid(True)

plt.xlabel("r /GeV")
plt.ylabel("sigma/sigma_0")
plt.title("Dipole Cross-Section Q2=50,  x=0.001")


plt.legend()
plt.savefig("/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots"+name+"/Model dependence 50 massive.png")
#plt.show()
plt.clf()
				
		
		
