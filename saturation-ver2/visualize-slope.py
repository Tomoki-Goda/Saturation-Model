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
    for j in ['2','4','6']:
        with open(args[0]+'/F2-slope-'+j+'.txt' ,"r") as fi:
            dpi=[]
            ri=[]
            for i in fi:
                data=i.strip().split("\t")
                dpi.append(float(data[1]))
                ri.append(float(data[0]))
        leg.append( ax1.plot(ri,dpi ,c='blue',ls="--"))
        ax1.text(ri[(len(ri)//5)*3],dpi[(len(ri)//5)*3],"x=$10^{-"+j+"}$" )
        with open(args[1]+'/F2-slope-'+j+'.txt',"r") as fi:
            dpi=[]
            ri=[]
            for i in fi:
                data=i.strip().split("\t")
                dpi.append(float(data[1]))
                ri.append(float(data[0]))
        leg.append(ax1.plot(ri,dpi ,c='red',ls="-"))
    
    ax1.legend([leg[0][0],leg[3][0]],['Without Sudakov','With Sudakov'])
    #ax1.set(title="",  ylabel="$\\frac{\\log\\partial F_2}{\\log\\partial x}$",    xlabel="$Q^2$",  xscale= 'log' ,   yscale='log' )
    ax1.set( xscale= 'log' ,   yscale='linear' )
    ax1.grid('true')
    #ax1.set_ylabel("$-\\frac{\\partial \\log F_2}{\\partial \\log x}$",rotation="vertical",loc='top')
    ax1.set_ylabel("$x-slope $",rotation="vertical",loc='top')
    ax1.set_xlabel("$Q^2$",rotation="horizontal",loc='right')
    fig1.set_figheight(5)
    fig1.set_figwidth(6)
    fig1.subplots_adjust(bottom=0.11, right=0.95, top=0.9, left=0.13)
    
    if saveflag:
        fig1.savefig(save1)
    else:
        plt.show()       
            
main()
    
