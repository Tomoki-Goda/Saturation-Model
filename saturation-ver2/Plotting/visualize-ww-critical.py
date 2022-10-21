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
    fig1,ax1=plt.subplots(1,2 ,sharey=True,layout='constrained')
    leg=[]
    #q2= ['100','650']
    #types= ['','tmd-','ww-']
    #names=['Curvature', 'Peak $\\alpha_s f$','Peak $\\Phi$']
    types= ['ww-']
    name=['GBW','BGK']

    for l in range(2):
        ax1[l].text(1.0e-6,0.4,name[l],fontsize=25)
        for j in range(len(types)):
            with open(args[0+2*l]+'/ww-critical-500.txt' ,"r") as fi:
                dpi=[]
                ri=[]
                for i in fi:
                    data=i.strip().split("\t")
                    dpi.append(float(data[1]))
                    ri.append(float(data[0]))
            leg.append( ax1[l].plot(ri,dpi ,c='blue',ls="--"))
            for k in ['5','50', '500']:    
                with open(args[1+2*l]+('/ww-critical-{0}.txt'.format(k)) ,"r") as fi:
                    dpi=[]
                    ri=[]
                    for i in fi:
                        data=i.strip().split("\t")
                        dpi.append(float(data[1]))
                        ri.append(float(data[0]))
                leg.append(ax1[l].plot(ri,dpi ,c='red',ls="-"))
                ax1[l].text(ri[0]*1.1,dpi[0],'$Q^2={0} \\mathrm{{GeV^2}}$'.format(k)) 
                ax1[l].set(xscale='log',yscale='log')
                ax1[l].grid(True)
    ax1[0].legend([leg[0][0],leg[1][0]],['Without Sudakov','With Sudakov'])
    #ax1.legend([leg[0][0],leg[3][0]],['Without Sudakov','With Sudakov'])
    #ax1.grid('true')
    ax1[0].set_ylabel("$Q^2_s\\;[\\mathrm{GeV^2}]$",rotation="horizontal",loc='top')
    ax1[1].set_xlabel("$x$",rotation="horizontal",loc='right')
    #fig1.subplots_adjust(bottom=0.1, right=0.95, top=0.95, left=0.1)
    fig1.set_figheight(4)
    fig1.set_figwidth(10)
    
    #fig1.subplots_adjust(bottom=0.1, right=0.95, top=0.95, left=0.1)
    
    if saveflag:
        plt.savefig(save1)
    else:
        plt.show()       
            
main()
    
