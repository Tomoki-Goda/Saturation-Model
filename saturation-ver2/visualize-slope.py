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
    fig1.get_layout_engine().set(wspace=1/10)
    leg=[]
    name=['GBW','BGK']
    for l in range(2):
        ax1[l].text(1.0e-2,0.35,name[l],fontsize=25)
        for j in ['2','4','6']:
            with open(args[0+2*l]+'/F2-slope-'+j+'.txt' ,"r") as fi:
                 dpi=[]
                 ri=[]
                 for i in fi:
                     data=i.strip().split("\t")
                     dpi.append(float(data[1]))
                     ri.append(float(data[0]))
            leg.append( ax1[l].plot(ri,dpi ,c='blue',ls="--"))
            ax1[l].text(ri[(len(ri)//5)*3],dpi[(len(ri)//5)*3],"x=$10^{-"+j+"}$" )
            with open(args[1+2*l]+'/F2-slope-'+j+'.txt',"r") as fi:
                dpi=[]
                ri=[]
                for i in fi:
                    data=i.strip().split("\t")
                    dpi.append(float(data[1]))
                    ri.append(float(data[0]))
            leg.append(ax1[l].plot(ri,dpi ,c='red',ls="-"))
        ax1[l].set( xscale= 'log' ,   yscale='linear' )
        ax1[l].grid('true')
        ax1[l].set_xlabel("$Q^2\\;[\\mathrm{{GeV^2}}]$",rotation="horizontal",loc='right')
    ax1[0].legend([leg[0][0],leg[3][0]],['Without Sudakov','With Sudakov'])
    #ax1.set_ylabel("$-\\frac{\\partial \\log F_2}{\\partial \\log x}$",rotation="vertical",loc='top')
    ax1[0].set_ylabel("$x-slope $",rotation="vertical",loc='top')
    fig1.set_figheight(4)
    fig1.set_figwidth(10)
    #fig1.subplots_adjust(bottom=0.11, right=0.95, top=0.9, left=0.13)
    
    if saveflag:
        fig1.savefig(save1)
    else:
        plt.show()       
            
main()
    
