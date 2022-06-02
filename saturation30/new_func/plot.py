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
	
################################## x dependence ###################################
dataarray=import_array("./Archive0106/GBSrfix/gbs500masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massless rfix log10x=-3")

dataarray=import_array("./Archive0106/GBSrfix/gbs500masslessLCB/pointsQ2500x5.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massless rfix log10x=-5")

dataarray=import_array("./Archive0106/GBSrfix/gbs10masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 10 massless rfix log10x=-3")

dataarray=import_array("./Archive0106/GBSrfix/gbs10masslessLCB/pointsQ2500x5.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 10 massless rfix log10x=-5")


plt.xscale('log')
plt.yscale('log')
plt.grid(True)

plt.xlabel("r /GeV")
plt.ylabel("sigma/sigma_0")
plt.title("GBS (r=0.5,m_l=0) Dipole Cross-Section Q^2=100GeV")


plt.legend()
plt.savefig('/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots/GBSrfix x dependence.png')
#plt.show()
plt.clf()
############################## Q dependence ####################################
dataarray=import_array("./Archive0106/GBSrfix/gbs500masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massless rfix Q2=500")

dataarray=import_array("./Archive0106/GBSrfix/gbs500masslessLCB/pointsQ2100x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massless rfix Q2=100")

dataarray=import_array("./Archive0106/GBSrfix/gbs10masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 10 massless rfix Q2=500")

dataarray=import_array("./Archive0106/GBSrfix/gbs10masslessLCB/pointsQ2100x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 10 massless rfix Q2=100")


plt.xscale('log')
plt.yscale('log')
plt.grid(True)

plt.xlabel("r /GeV")
plt.ylabel("sigma/sigma_0")
plt.title("GBS (r=0.5,m_l=0) Dipole Cross-Section x=0.001")


plt.legend()
plt.savefig('/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots/GBSrfix Q2 dependence.png')
#plt.show()
plt.clf()	

############################## Model dependence ####################################
dataarray=import_array("./Archive0106/GBS/gbs500masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massless")


dataarray=import_array("./Archive0106/GBSrfix/gbs500masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massless rmax=0.5")

dataarray=import_array("./Archive0106/GBSPert/gbs500masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 500 massless")

dataarray=import_array("./Archive0106/GBSPertrfix/gbs500masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 500 massless rmax=0.5")

dataarray=import_array("./GBWIntegratedHighAccu/500masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 500 massless")

plt.xscale('log')
plt.yscale('log')
plt.grid(True)

plt.xlabel("r /GeV")
plt.ylabel("sigma/sigma_0")
plt.title("Dipole Cross-Section Q2=500,  x=0.001")


plt.legend()
plt.savefig('/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots/Model dependence 500 massless.png')
#plt.show()
plt.clf()

#########################################massive 500 #######################################################
dataarray=import_array("./Archive0106/GBS/gbs500massiveLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massive")

dataarray=import_array("./Archive0106/GBSrfix/gbs500massiveLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 500 massive rmax=0.5")

dataarray=import_array("./Archive0106/GBSPert/gbs500massiveLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 500 massive")

dataarray=import_array("./Archive0106/GBSPertrfix/gbs500massiveLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 500 massive rmax=0.5")

dataarray=import_array("./GBWIntegratedHighAccu/500massiveLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 500 massive")


plt.xscale('log')
plt.yscale('log')
plt.grid(True)

plt.xlabel("r /GeV")
plt.ylabel("sigma/sigma_0")
plt.title("Dipole Cross-Section Q2=500,  x=0.001")


plt.legend()
plt.savefig('/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots/Model dependence 500 massive.png')
#plt.show()
plt.clf()

############################################ massless 100 #################################################################	
dataarray=import_array("./Archive0106/GBS/gbs100masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 100 massless")

dataarray=import_array("./Archive0106/GBSrfix/gbs100masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 100 massless rmax=0.5")

dataarray=import_array("./Archive0106/GBSPert/gbs100masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 100 massless")

dataarray=import_array("./Archive0106/GBSPertrfix/gbs100masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 100 massless rmax=0.5")

dataarray=import_array("./GBWIntegratedHighAccu/100masslessLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 100 massless")

plt.xscale('log')
plt.yscale('log')
plt.grid(True)

plt.xlabel("r /GeV")
plt.ylabel("sigma/sigma_0")
plt.title("Dipole Cross-Section Q2=500,  x=0.001")


plt.legend()
plt.savefig('/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots/Model dependence 100 massless.png')
#plt.show()
plt.clf()

############################################ massive 100 ####################################################################
dataarray=import_array("./Archive0106/GBS/gbs100massiveLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 100 massive")

dataarray=import_array("./Archive0106/GBSrfix/gbs100massiveLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs 100 massive rmax=0.5")

dataarray=import_array("./Archive0106/GBSPert/gbs100massiveLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 100 massive")

dataarray=import_array("./Archive0106/GBSPertrfix/gbs100massiveLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbs pert 100 massive rmax=0.5")

dataarray=import_array("./GBWIntegratedHighAccu/100massiveLCB/pointsQ2500x3.txt")
plt.plot(dataarray[0],dataarray[1],label="gbw 100 massive")


plt.xscale('log')
plt.yscale('log')
plt.grid(True)

plt.xlabel("r /GeV")
plt.ylabel("sigma/sigma_0")
plt.title("Dipole Cross-Section Q2=500,  x=0.001")


plt.legend()
plt.savefig('/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes/NewPlots/Model dependence 100 massive.png')
#plt.show()
plt.clf()
				
		
		
