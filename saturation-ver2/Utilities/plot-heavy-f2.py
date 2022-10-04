#! /usr/bin/env python3
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import math
import sys 
import getopt

#from matplotlib.ticker import LogLocator
import matplotlib.ticker as tk
#fig,ax=plt.subplots()
#mpl.rcParams["xaxis.labellocation"]="right"
#mpl.rcParams["yaxis.labellocation"]="top"
 
#data_array=[]
x=[]
q2=[]
f2=[]
x2=[]
q22=[]
f22=[]

xd=[]
q2d=[]
f2d=[]
f2e=[]
f2c=[]
counter=0
def shift_col(val, low, high):
    return( (val-low)/(high-low))
    
    
def sortingkey(datalist):
    return(float(datalist[0]))
    
    
def main():
    #data=[]
    #res=[]
    saveflag=False
    limit_plot=False
    col1=3
    name="F_2"
    argv=sys.argv[1:]
    try:
        opts , args =getopt.getopt(argv,"s:lcb" )
        labels=args;
        plotstyle=["-" for i in range(len(args))]
    except:
        print("option eror");
        
    for opt,arg in opts:
        if opt in ["-s","--save"]:
            save1=arg+"1"
            save2=arg+"2"
            save3=arg+"hist"
            save4=arg+"hist-diff"
            saveflag=True
            print("SAVE")
        if opt in ['-l','--limit']:
            limit_plot=True
        if opt in ['-c','--charm']:
            col1=3
            name="F_2^c"
        if opt in ['-b','--bottom']:
            col1=2
            name="F_2^b"
        else:
            print("Unknown") 
        
        #opts , args =getopt.getopt(argv,"")
    #for i in args:
    
    with open(args[0]+'.txt',"r") as fi:
        for i in fi:
            data=i.strip().split("\t")
            #arr.append(data)
            x.append(data[1])
            f2.append(data[2])
            q2.append(data[0])
    #print(len(q2))
    with open(args[1]+'.txt',"r") as fi:
        for i in fi:
            data=i.strip().split("\t")
            #arr.append(data)
            x2.append(data[1])
            f22.append(data[2])
            q22.append(data[0])
            
    with open(args[0]+'-data.txt',"r") as fi:
        arr=fi.readlines();
    print(len(arr))
    
    arr=[i.strip().split("\t") for i in arr]
    
    arr=np.transpose(arr)
    f2d=arr[2]
    f2e=arr[3]
    xd=arr[1]
    q2d=arr[0]
    with open(args[1]+'-data.txt',"r") as fi:
        arr=fi.readlines();
    print(len(arr))
    
    arr=[i.strip().split("\t") for i in arr]
    
    arr=np.transpose(arr)
    f2e2=arr[3]
    f2cd=arr[2]
    xd2=arr[1]
    q2d2=arr[0]
    
    
    q2set=sorted( list(set(q2)) , key=float)
    q2len=len(q2set)
    print(q2len, " frames")
    #col1=3
    row1=(q2len)//col1 
    if((q2len%col1)!=0):
         row1+=1
    
    print(row1,"*", col1)
     
    fig1,ax1=plt.subplots( nrows=row1, ncols=col1,sharex="col",sharey="row",constrained_layout=True)

    fig1.set_constrained_layout_pads(hspace=0,wspace=0)
    ######################################################################################################
    legend=[]
    for i in range(q2len):
        pos=[i//col1,i%col1]
        q2val=q2set[i];
        #print(q2val);
        arr=[ ]
        arr2=[ ]
        #yarr=[]
        for j in range(len(q2)):
            if(q2[j]==q2set[i]):
                #print(q2set[i])
                arr.append([ x[j],f2[j] ])
                arr2.append([ x2[j],f22[j] ])
        arr=np.transpose(sorted(arr,key=sortingkey))
        #print(arr)
        arr2=np.transpose(sorted(arr2,key=sortingkey))
        #print( arr)
        xarr=[float(k) for k in arr[0]]
        yarr=[float(k) for k in arr[1]]
        xarr2=[float(k) for k in arr2[0]]
        yarr2=[float(k) for k in arr2[1]]
        #print( xarr, yarr)

        legend.append(ax1[pos[0]][pos[1]].plot(xarr,yarr ,c='blue',ls="--"))
        legend.append(ax1[pos[0]][pos[1]].plot(xarr2,yarr2 ,c='red'))    
        #########################################################################            
        for j in range(len(q2d)):
            if(q2d[j]==q2set[i]):
                ax1[pos[0]][pos[1]].errorbar(float(xd[j]),float(f2d[j]),yerr=float(f2e[j]),c='green')
                ax1[pos[0]][pos[1]].scatter(float(xd[j]),float(f2d[j]),s=12,marker="x",c='green')
        #########################################################################        
        ax1[pos[0]][pos[1]].set(  xscale="log" ,   yscale='linear')
        ax1[pos[0]][pos[1]].text(xarr[0]*0.8, yarr[0],"$Q^2={val}\\;\\mathrm{{GeV^2}}$".format(val=float(q2val)) , fontsize=8)
        ax1[pos[0]][pos[1]].yaxis.set_major_locator(plt.LogLocator(base=10))
        ax1[pos[0]][pos[1]].yaxis.set_minor_locator(plt.NullLocator())
        ax1[pos[0]][pos[1]].xaxis.set_major_locator(plt.FixedLocator([1.0e-3,1.0e-5] ))
        ax1[pos[0]][pos[1]].xaxis.set_minor_locator(plt.NullLocator())
        ax1[pos[0]][pos[1]].margins(x=0)
####################################################################################
    ax1[0][0].legend([legend[0][0],legend[1][0]],["Without Sudakov", "With Sudakov"])
    fig1.set_figheight(2*row1+1)
    fig1.set_figwidth(2*col1+1)
    fig1.supylabel("$"+name+"$",rotation="horizontal",x=0.02,y=0.95,fontsize=13)
    fig1.supxlabel("$x$",rotation="horizontal",y=0.015,x=0.99,fontsize=13)
    if saveflag:
        fig1.savefig(save1)
    else:
    	plt.show()
        
    return(0)
main()




















