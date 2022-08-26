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
    fig1,ax1=plt.subplots(1,2 ,sharey=True,constrained_layout=True)
    leg=[]
    #q2= ['100','650']
    types= ['','tmd-']
    names=['Curvature', 'Peak']
    for j in range(len(types)):
        with open(args[0]+'/'+types[j]+'critical-650.txt' ,"r") as fi:
            dpi=[]
            ri=[]
            for i in fi:
                data=i.strip().split("\t")
                dpi.append(float(data[1]))
                ri.append(float(data[0]))
        leg.append( ax1[j].plot(ri,dpi ,c='blue',ls="--"))
        ax1[j].set(xscale="log" ,   yscale='linear' )
        ax1[j].grid('true')
        
        
        with open(args[1]+'/'+types[j]+'critical-650.txt' ,"r") as fi:
            dpi=[]
            ri=[]
            for i in fi:
                data=i.strip().split("\t")
                dpi.append(float(data[1]))
                ri.append(float(data[0]))
        leg.append(ax1[j].plot(ri,dpi ,c='red',ls="-"))
        ax1[j].set(xscale="log" ,   yscale='log' )
        ax1[j].grid('true')
        ax1[j].text(ri[len(dpi)//2]*2,dpi[len(dpi)//2]*1.5,names[j])
        
    ax1[1].legend([leg[0][0],leg[1][0]],['Without Sudakov','With Sudakov'])
    #ax1.legend([leg[0][0],leg[3][0]],['Without Sudakov','With Sudakov'])
    #ax1.set( xscale= 'log' ,   yscale='log' )
    #ax1.grid('true')
    ax1[0].set_ylabel("$Q^2_s$",rotation="horizontal",loc='top')
    ax1[1].set_xlabel("$x$",rotation="horizontal",loc='right')
    #fig1.subplots_adjust(bottom=0.1, right=0.95, top=0.95, left=0.1)
    fig1.set_figheight(4)
    fig1.set_figwidth(8)
    if saveflag:
        fig1.savefig(save1)
    else:
        plt.show()       
            
main()
    
