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
    fig1,ax1=plt.subplots( )
    leg=[]
    #for k in ['100','650']:
    for j in ['1']:
        with open(args[0]+'/gluon-x-650-'+j+'.txt' ,"r") as fi:
            dpi=[]
            ri=[]
            for i in fi:
                data=i.strip().split("\t")
                dpi.append(float(data[1]))
                ri.append(float(data[0]))
        leg.append( ax1.plot(ri,dpi ,c='blue',ls="--"))
                    
        with open(args[1]+'/gluon-x-100-'+j+'.txt',"r") as fi:
            dpi=[]
            ri=[]
            for i in fi:
                data=i.strip().split("\t")
                dpi.append(float(data[1]))
                ri.append(float(data[0]))
        leg.append(ax1.plot(ri,dpi ,c='red',ls="-"))
        
        
        with open(args[1]+'/gluon-x-650-'+j+'.txt' ,"r") as fi:
            dpi=[]
            ri=[]
            for i in fi:
                data=i.strip().split("\t")
                dpi.append(float(data[1]))
                ri.append(float(data[0]))
        leg.append( ax1.plot(ri,dpi ,c='magenta',ls="-"))
        
    
    
    ax1.legend([leg[0][0],leg[1][0],leg[2][0]],['Without Sudakov','With Sudakov $Q^2=100\\;\\mathrm{GeV}^2$','With Sudakov $Q^2=650\\;\\mathrm{GeV}^2$'])
    #ax1.set(title="",  ylabel="$\\alpha_s f(x k^2)$",    xlabel="$k^2$",  xscale= 'log' ,   yscale='linear' )
    ax1.set( xscale= 'log' ,   yscale='linear' )
    ax1.grid('true')
    ax1.set_ylabel('$\\alpha_s f(x k^2)$',rotation="vertical",loc='top')
    ax1.set_xlabel("$x$",rotation="horizontal",loc='right')
    fig1.set_figheight(5)
    fig1.set_figwidth(6)
    fig1.subplots_adjust(bottom=0.1, right=0.95, top=0.95, left=0.1)
    if saveflag:
        fig1.savefig(save1)
    else:
        plt.show()              
            
main()
