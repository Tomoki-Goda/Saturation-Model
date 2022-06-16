#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
#import math
savedir="/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes"

def sudakov_np(r, Q2):
	g1=0.01
	g2=0.2
	Q_02=1
	r_max=0.5
	val= g1/2*pow(r,2)+np.log(1+pow(r/r_max,2)) *np.log(Q2/Q_02) *g2/4 
	return(val)
	
def sudakov_p(r, Q2):
	r_max=0.5
	m2=1.26*(pow(r,-2)+pow(r_max,-2))
	if( m2 > Q2):
		return(0)
	
	logQ=np.log(Q2/0.09)
	logm=np.log(m2/0.09)
	
	val=logQ*np.log(logQ/logm)-logQ+logm
	return(val)
	
	
		

x=[ 0.1973*pow(10,-2+3*(i/100))  for i in range(100)]
#print(x)

#ynp=[sudakov_np(i,Q2) for i in x]

#yp=[sudakov_p(i,Q2) for i in x]

#ypnp=[sudakov_p(i,Q2) +sudakov_np(i,Q2) for i in x]
Q2=5
ynp=[np.exp(-sudakov_np(i,Q2)) for i in x]
yp=[np.exp(-sudakov_p(i,Q2)) for i in x]
ypnp=[np.exp(-(sudakov_p(i,Q2) +sudakov_np(i,Q2))) for i in x]

plt.plot(x,yp,label="Q2=5 Gev pert",ls="-")
plt.plot(x,ynp,label="Q2=5 non-pert",ls="-.")
plt.plot(x,ypnp,label="Q2=5 total",ls=":")

Q2=500
ynp=[np.exp(-sudakov_np(i,Q2)) for i in x]
yp=[np.exp(-sudakov_p(i,Q2)) for i in x]
ypnp=[np.exp(-(sudakov_p(i,Q2) +sudakov_np(i,Q2))) for i in x]

plt.plot(x,yp,label="Q2=500 Gev pert",ls="-")
plt.plot(x,ynp,label="Q2=500 Gev non-pert",ls="-.")
plt.plot(x,ypnp,label="Q2=500 Gev total",ls=":")


plt.title("Exp(-Sudakov)")

plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.xlabel("r (integration variable!)")
plt.ylabel("e^-S")

plt.savefig(savedir+"/Sudakov.png")
plt.show()
plt.clf()


