#! /usr/bin/env python3 
import numpy as np
import matplotlib.pyplot as plt

savedir="/media/tomoki/TOMOKI-USB/Saturation Model/Saturation Notes"

def QuickPlot(func, llim, ulim,plt,plotlabel="none",linestyle="-"):
	plotnum=1000
	interval=(ulim-llim)/plotnum
	arrx=[llim+i*interval for i in range(plotnum)]
	arry=[func(i) for i in arrx]
	plt.plot(arrx,arry,label=plotlabel,ls=linestyle)
	plt.grid(True)
	#plt.xscale("log")
	#plt.yscale("log")
	#plt.show()
	
	
def sudakov_np(r, Q2):
	g1=0.013
	g2=0.2
	Q_0=1
	r_max=0.5
	val= g1/2*pow(r,2)+np.log(1+pow(r/r_max,2)) *np.log(Q2/Q_0) *g2/4 
	return(val)
	
def sudakov_p(r, Q2):
	r_max=0.5
	m2=1.26*(pow(r,-2)+pow(r_max,-2))
#	if( m2 > pow(Q,2)):
#		return(0)
	
	logQ=np.log(Q2/0.09)
	logm=np.log(m2/0.09)
	
	val=logQ*np.log(logQ/logm)-logQ+logm
	return(val)
	
	

	
def integrand(r,R):
	#R=5
	l=0.28
	x0=0.0001
	x=0.001
	Qs2=pow(x0/x,l)
	s0=23.0
	
	val=r*np.log(R/r)*s0 * Qs2*(1-Qs2*pow(r/2,2))*np.exp(-pow(r/2,2)*Qs2) 
	return(val)
	

def int5(r):
	return(integrand(r,5))

def integrandSnp5(r):
	Q=10
	val=int5(r)*np.exp(-sudakov_p(r,Q)-sudakov_np(r,Q))
	return(val)
	
def integrandS5(r):
	Q=10
	val=int5(r)*np.exp(-sudakov_p(r,Q))
	return(val)
	
def int1(r):
	return(integrand(r,1))
	
def integrandSnp1(r):
	Q=10
	val=int1(r)*np.exp(-sudakov_p(r,Q)-sudakov_np(r,Q))
	return(val)
	
def integrandS1(r):
	Q=10
	val=int1(r)*np.exp(-sudakov_p(r,Q))
	return(val)

def int01(r):
	return(integrand(r,0.1))
	
def integrandSnp01(r):
	Q=10
	val=int01(r)*np.exp(-sudakov_p(r,Q)-sudakov_np(r,Q))
	return(val)
	
def integrandS01(r):
	Q=10
	val=int01(r)*np.exp(-sudakov_p(r,Q))
	return(val)

QuickPlot(integrandSnp01,0.0001,0.1,plt,"SudNP","-")
QuickPlot(integrandS01,0.0001,0.1,plt,"SudPert","-.")
QuickPlot(int01,0.0001,0.1,plt,"GBW",":")
QuickPlot(integrandSnp1,0.0001,1,plt,"SudNP","-")
QuickPlot(integrandS1,0.0001,1,plt,"SudPert","-.")
QuickPlot(int1,0.0001,1,plt,"GBW",":")
QuickPlot(integrandSnp5,0.0001,5,plt,"SudNP","-")
QuickPlot(integrandS5,0.0001,5,plt,"SudPert","-.")
QuickPlot(int5,0.0001,5,plt,"GBW",":")

#plt.xscale("log")

plt.xlabel("r GeV")
plt.ylabel("integrand")
plt.title("integrand at different upper limit r Q2=10,x=0.0001")
plt.legend()
plt.savefig(savedir+"/IntegrandsQ210.png")
plt.show()






















