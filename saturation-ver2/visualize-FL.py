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
    fig1,ax1=plt.subplots(1,2,constrained_layout=True,sharex=True, sharey=True )
    fig1.get_layout_engine().set(wspace=1/10)
    leg=[]
    name=['GBW','BGK']
    for l in range(2):
        ax1[l].text(0.75,3.5e-1,name[l],fontsize=25)
        #for w in ['100', '200', '300']:   
        for w in ['200','data']:   
            with open(args[0+2*l]+'/FL-'+w+'.txt' ,"r") as fi:
                 dpi=[]
                 ri=[]
                 for i in fi:
                     data=i.strip().split("\t")
                     dpi.append(float(data[1]))
                     ri.append(float(data[0]))
            if w=='data':
                leg.append(ax1[l].step(ri,dpi,c="blue", ls="--",where='mid')) 
            else:
                leg.append( ax1[l].plot(ri,dpi ,c='blue',ls="--"))
            #ax1[l].text(float(data[0])/2,float(data[1])+0.01,"W="+w+'GeV')
             
            with open(args[1+2*l]+'/FL-'+w+'.txt',"r") as fi:
                dpi=[]
                ri=[]
                for i in fi:
                    data=i.strip().split("\t")
                    dpi.append(float(data[1]))
                    ri.append(float(data[0]))
            if w=='data':
                leg.append(ax1[l].step(ri,dpi,c='red',ls="-",where="mid"))
            else:
                leg.append(ax1[l].plot(ri,dpi ,c='red',ls="-"))
###########################################            
        with open('./data/H1FL.txt','r') as fi:
            expdata=fi.readlines()
            expdata=[i.strip().split() for i in expdata]
            #expdata=np.transpose(expdata)
            w=[np.sqrt(float(i[0])*(1-float(i[1]))/float(i[1])) for i in expdata]
            #print(w)
            expdata=[ [float(i[0]) for i in expdata], [float(i[2]) for i in expdata], [float(i[6]) for i in expdata]]
            
        #print(expdata)
        ax1[l].scatter(expdata[0],expdata[1],c='g')
        for i in range(len(expdata[0])):
            #print(expdata[0][i],expdata[1][i],expdata[2][i])
            ax1[l].errorbar(expdata[0][i],expdata[1][i],yerr=expdata[2][i],c='g')
            if i%2==0:
                ax1[l].text(expdata[0][i],-0.025,'{0:.0f}'.format(w[i]),rotation=70)
            else:
                ax1[l].text(expdata[0][i],-0.1,'{0:.0f}'.format(w[i]),rotation=70)
                
        ax1[l].set( xscale= 'log' ,   yscale='linear' )
        ax1[l].grid('true')
        ax1[l].text(0.75,-0.05,'$W$=',fontsize=8)
        ax1[l].set_xlabel("$Q^2 \\;[\\mathrm{GeV^2}]$",rotation="horizontal",loc='right')
#################################################
    ax1[0].legend([leg[0][0],leg[1][0]],['Without Sudakov','With Sudakov'])
    ax1[0].set_ylabel("$F_L$",rotation="vertical",loc='top')
    #plt.xlim([1,1.0e+3])
    plt.ylim([-0.2,0.5])
    fig1.set_figheight(4)
    fig1.set_figwidth(10.5)
    #fig1.subplots_adjust(bottom=0.11, right=0.95, top=0.95, left=0.11)
    if saveflag:
        fig1.savefig(save1)
    else:
        plt.show()             
            
main()
