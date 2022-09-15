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
    leg=[]
    name=['GBW','BGK']
    for l in range(2):
        ax1[l].text(1.0e-2,3.5,name[l],fontsize=25);
        for j in ['4']:
            with open(args[0+2*l]+'/gluon-500-'+j+'.txt' ,"r") as fi:
                dpi=[]
                ri=[]
                for i in fi:
                    data=i.strip().split("\t")
                    dpi.append(float(data[1]))
                    ri.append(float(data[0]))
            leg.append( ax1[l].plot(ri,dpi ,c='blue',ls="--"))
            q2=['5','50','500']
            #sty=['-','_.',':']
            for k in range(len(q2)):
                #print(k)
                with open(args[1+2*l]+('/gluon-{0}-'.format(q2[k]))+j+'.txt',"r") as fi:
                    dpi=[]
                    ri=[]
                    for i in fi:
                        data=i.strip().split("\t")
                        dpi.append(float(data[1]))
                        ri.append(float(data[0]))
                    pos=int(len(dpi)/2)
                    ax1[l].text(ri[pos]*3,dpi[pos], '$Q^2={0}\\;\\mathrm{{ GeV^2}}$'.format(q2[k]) )
                leg.append(ax1[l].plot(ri,dpi ,c='red',ls="-"))

        ax1[l].set( xscale= 'log' ,   yscale='linear' )
        ax1[l].grid('true')
#    ax1[0].legend([leg[0][0],leg[1][0],leg[2][0]],['Without Sudakov','With Sudakov $Q^2=100\\;\\mathrm{GeV}^2$','With Sudakov $Q^2=650\\;\\mathrm{GeV}^2$'])
    ax1[0].legend([leg[0][0],leg[1][0]],['Without Sudakov','With Sudakov'])
    #ax1.set(title="",  ylabel="$\\alpha_s f(x k^2)$",    xlabel="$k^2$",  xscale= 'log' ,   yscale='linear' )
    ax1[0].set_ylabel('$\\alpha_s \\mathcal{F}(x, k^2)$',rotation="vertical",loc='top')
    ax1[1].set_xlabel("$k^2\\;(\\mathrm{GeV}^2) $",rotation="horizontal",loc='right')
    fig1.set_figheight(4)
    fig1.set_figwidth(10)
    #fig1.subplots_adjust(bottom=0.1, right=0.95, top=0.95, left=0.1)
    if saveflag:
        fig1.savefig(save1)
    else:
        plt.show()              
            
main()
