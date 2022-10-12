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
    fig1,ax1=plt.subplots(1,2,constrained_layout=True,sharey=True,sharex=True )
    fig1.get_layout_engine().set( wspace=1/10)
    leg=[]
    name=['GBW','BGK']
    for l in range(2): 
        ax1[l].text(1,1.0e-5,name[l],fontsize=25)
        for j in ['2', '4', '6']:
            with open(args[0+2*l]+'/dipole-100-'+j+'.txt' ,"r") as fi:
                dpi=[]
                ri=[]
                for i in fi:
                    data=i.strip().split("\t")
                    dpi.append(float(data[1]))
                    ri.append(float(data[0]))
                
            leg.append( ax1[l].plot(ri,dpi ,c='blue',ls="--"))
            ax1[l].text(ri[0]/1.3,5*dpi[0],"x=$10^{-"+j+"}$")
            with open(args[1+2*l]+'/dipole-100-'+j+'.txt',"r") as fi:
                dpi=[]
                ri=[]
                for i in fi:
                    data=i.strip().split("\t")
                    dpi.append(float(data[1]))
                    ri.append(float(data[0]))
            leg.append(ax1[l].plot(ri,dpi ,c='red',ls="-"))
        ax1[l].set( xscale= 'log' ,   yscale='log' )
        ax1[l].grid('true',which='major')
        ax1[l].set_xlabel("$r\\;[\\mathrm{GeV}^{-1}]$",rotation="horizontal",loc='right')
    ax1[0].legend([leg[0][0],leg[3][0]],['Without Sudakov','With Sudakov'])
    ax1[0].set_ylabel("$\\sigma/\\sigma_0$",rotation="vertical",loc='top')
    fig1.set_figheight(4)
    fig1.set_figwidth(10.5)
    #fig1.subplots_adjust(bottom=0.11, right=0.95, top=0.95, left=0.11)
    if saveflag:
        fig1.savefig(save1)
    else:
        plt.show()      
    
    
            
main()
    



