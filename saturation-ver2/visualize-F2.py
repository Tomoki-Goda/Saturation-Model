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
    fig1,ax1=plt.subplots(1,2,constrained_layout=True,sharex=True,sharey=True )
    leg=[]
    name=['GBW', 'BGK']
    for l in range(2): 
        ax1[l].text(1.0e+1,1.0e-2,name[l],fontsize=25);
        for j in ['2','4','6']:
            with open(args[0+2*l]+'/F2-'+j+'.txt' ,"r") as fi:
                dpi=[]
                ri=[]
                for i in fi:
                    data=i.strip().split("\t")
                    dpi.append(float(data[1]))
                    ri.append(float(data[0]))
            leg.append( ax1[l].plot(ri,dpi ,c='blue',ls="--"))
            ax1[l].text(float(data[0])/2,float(data[1]),"x=$10^{-"+j+"}$")
            
            with open(args[1+2*l]+'/F2-'+j+'.txt',"r") as fi:
                dpi=[]
                ri=[]
                for i in fi:
                    data=i.strip().split("\t")
                    dpi.append(float(data[1]))
                    ri.append(float(data[0]))
            leg.append(ax1[l].plot(ri,dpi ,c='red',ls="-"))
        ax1[l].set( xscale= 'log' ,   yscale='log' )
        ax1[l].grid('true')
    
    ax1[0].legend([leg[0][0],leg[3][0]],['Without Sudakov','With Sudakov'])
    ax1[0].set_ylabel("$F_2$",rotation="vertical",loc='top',fontsize=13)
    ax1[1].set_xlabel("$Q^2 \\;(\\mathrm{GeV^2})$",rotation="horizontal",loc='right',fontsize=13)
    fig1.set_figheight(5)
    fig1.set_figwidth(10)
    #fig1.subplots_adjust(bottom=0.11, right=0.95, top=0.95, left=0.11)
    if saveflag:
        fig1.savefig(save1)
    else:
        plt.show()             
            
main()
