#! /usr/bin/env python3
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import math
import sys 
import getopt

#from matplotlib.ticker import LogLocator
import matplotlib.ticker as tk

def main():
    #data=[]
    saveflag=False
    argv=sys.argv[1:]
    try:
        opts , args =getopt.getopt(argv,"s:" )
        labels=args;
        plotstyle=["-" for i in range(len(args))]
    except:
        print("option eror");
        
    for opt,arg in opts:
        if opt in ["-s","--save"]:
            save1=arg
            saveflag=True
            print("SAVE")
        else:
            print("Unknown") 
        
    
    dp=[]
    dpi=[]
    r=[]
    ri=[]
    fig1,ax1=plt.subplots(1,1 ,sharey=True,layout='constrained')
    leg=[]
    #q2= ['100','650']
    #types= ['','tmd-','ww-']
    #names=['Curvature', 'Peak $\\alpha_s f$','Peak $\\Phi$']
    types= ['ww-']
    names=['Peak $\\Phi$']
    for j in range(len(types)):
        with open(args[0]+'/ww-critical-650.txt' ,"r") as fi:
            dpi=[]
            ri=[]
            for i in fi:
                data=i.strip().split("\t")
                dpi.append(float(data[1]))
                ri.append(float(data[0]))
        leg.append( ax1.plot(ri,dpi ,c='blue',ls="--"))
 
        with open(args[1]+'/ww-critical-5.txt' ,"r") as fi:
            dpi=[]
            ri=[]
            for i in fi:
                data=i.strip().split("\t")
                dpi.append(float(data[1]))
                ri.append(float(data[0]))
        leg.append(ax1.plot(ri,dpi ,c='red',ls="-"))
        ax1.text(ri[0]*1.1,dpi[0],'$Q^2=5 \mathrm{GeV^2}$')
       
        with open(args[1]+'/ww-critical-100.txt' ,"r") as fi:
            dpi=[]
            ri=[]
            for i in fi:
                data=i.strip().split("\t")
                dpi.append(float(data[1]))
                ri.append(float(data[0]))
        leg.append(ax1.plot(ri,dpi ,c='red',ls="-"))
        ax1.text(ri[0]*1.1,dpi[0],'$Q^2=100 \mathrm{GeV^2}$')
 
        
        with open(args[1]+'/ww-critical-650.txt' ,"r") as fi:
            dpi=[]
            ri=[]
            for i in fi:
                data=i.strip().split("\t")
                dpi.append(float(data[1]))
                ri.append(float(data[0]))
        leg.append(ax1.plot(ri,dpi ,c='red',ls="-"))
        ax1.set(xscale="log" ,   yscale='log' )
        ax1.grid('true')

        ax1.text(ri[0]*1.1,dpi[0],'$Q^2=650 \mathrm{GeV^2}$')
        #ax1.text(ri[len(dpi)//2]*2,dpi[len(dpi)//2]*1.5,names[j])
        
    ax1.legend([leg[0][0],leg[1][0]],['Without Sudakov','With Sudakov'])
    #ax1.legend([leg[0][0],leg[3][0]],['Without Sudakov','With Sudakov'])
    #ax1.set( xscale= 'log' ,   yscale='log' )
    #ax1.grid('true')
    ax1.set_ylabel("$Q^2_s$",rotation="horizontal",loc='top')
    ax1.set_xlabel("$x$",rotation="horizontal",loc='right')
    #fig1.subplots_adjust(bottom=0.1, right=0.95, top=0.95, left=0.1)
    fig1.set_figheight(5)
    fig1.set_figwidth(6)
    
    fig1.subplots_adjust(bottom=0.1, right=0.95, top=0.95, left=0.1)
    
    if saveflag:
        plt.savefig(save1)
    else:
        plt.show()       
            
main()
    
