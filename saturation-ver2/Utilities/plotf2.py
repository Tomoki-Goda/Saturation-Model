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
counter=0;
def shift_col(val, low, high):
    return( (val-low)/(high-low))
    
    
def sortingkey(datalist):
    return(float(datalist[0]))
    
    
def main():
    #data=[]
    #res=[]
    saveflag=False
    limit_plot=False
    argv=sys.argv[1:]
    try:
        opts , args =getopt.getopt(argv,"s:l" )
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
        else:
            print("Unknown") 
        
        #opts , args =getopt.getopt(argv,"")
    #for i in args:
    
    with open(args[0]+'/fcn.txt',"r") as fi:
        for i in fi:
            data=i.strip().split("\t")
            #arr.append(data)
            x.append(data[0])
            f2.append(data[1])
            q2.append(data[2])
    #print(len(q2))
    with open(args[1]+'/fcn.txt',"r") as fi:
        for i in fi:
            data=i.strip().split("\t")
            #arr.append(data)
            x2.append(data[0])
            f22.append(data[1])
            q22.append(data[2])
            
    with open(args[0]+'/data.txt',"r") as fi:
        arr=fi.readlines();
    print(len(arr))
    
    arr=[i.strip().split("\t") for i in arr]
    
    arr=np.transpose(arr)
    f2c=arr[0]
    f2d=arr[1]
    f2e=arr[2]
    xd=arr[3]
    q2d=arr[4]
    with open(args[1]+'/data.txt',"r") as fi:
        arr=fi.readlines();
    print(len(arr))
    
    arr=[i.strip().split("\t") for i in arr]
    
    arr=np.transpose(arr)
    f2c2=arr[0]
    f2d2=arr[1]
    f2e2=arr[2]
    xd2=arr[3]
    q2d2=arr[4]
    
    chi1=[]
    chi2=[]
    
    if limit_plot:
        #print(q2)
        q2set=["1.10000e-01","5.00000e-01", "6.50000e+00","1.80000e+01","4.50000e+01","5.00000e+02"]
        print(q2set)
        q2len=2*len(q2set)
        extra=0
        row1=2
        col1=3
        fig1,ax1=plt.subplots( nrows=row1, ncols=col1,sharex="col",sharey="row",constrained_layout=True)
        fig1.set_constrained_layout_pads(hspace=0,wspace=0)
    else:
        q2set=sorted( list(set(q2)) , key=float)
        q2len=len(q2set)
        col1=4
        row1=(q2len//2)//col1 
        extra=0
        if (((q2len//2))%col1 )!= 0 :
            row1+=1
        extra=col1*row1-(q2len//2)
        
        fig1,ax1=plt.subplots( nrows=row1, ncols=col1,sharex="col",sharey="row",constrained_layout=True)
        col2=4
        row2=(q2len-(q2len//2+extra))//col2
        if ( (q2len-(q2len//2+extra))%col2 )!= 0 :
            row2+=1
        
        fig2,ax2=plt.subplots( nrows=row2, ncols=col2,sharex=True,sharey=True,constrained_layout=True)

        fig1.set_constrained_layout_pads(hspace=0,wspace=0)
        fig2.set_constrained_layout_pads(hspace=0,wspace=0)
    ######################################################################################################
    for i in range(q2len//2 +extra):
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
        ax1[pos[0]][pos[1]].plot(xarr,yarr ,c='blue',ls="--")
        ax1[pos[0]][pos[1]].plot(xarr2,yarr2 ,c='red')    
        #########################################################################            
        chi1.append([])
        chi2.append([])
        for j in range(len(q2d)):
            if(q2d[j]==q2set[i]):
                ax1[pos[0]][pos[1]].errorbar(float(xd[j]),float(f2d[j]),yerr=float(f2e[j]),c='green')
                ax1[pos[0]][pos[1]].scatter(float(xd[j]),float(f2d[j]),s=12,marker="x",c='green')
                chi1[i].append(pow( (float(f2d[j])-float(f2c[j]))/float(f2e[j]) ,2))
            if(q2d2[j]==q2set[i]):
                chi2[i].append(pow( (float(f2d2[j])-float(f2c2[j]))/float(f2e2[j]) ,2))
        #########################################################################        
        #chi+=pow((float(f2d[j])-float(f2c[j]))/float(f2e[j]),2)
        #count+=1
        ax1[pos[0]][pos[1]].set(  xscale="log" ,   yscale='linear')
        ax1[pos[0]][pos[1]].text(xarr[0]*1.2, yarr[0],"$Q^2={val}\\;\\mathrm{{GeV^2}}$".format(val=float(q2val)),fontsize=7 )
        #ax1[pos[0]][pos[1]].yaxis.set_major_locator(plt.FixedLocator([0.01,0.05,0.1,0.5,1,2.5]))
        ax1[pos[0]][pos[1]].yaxis.set_major_locator(plt.LogLocator(base=10))
        ax1[pos[0]][pos[1]].yaxis.set_minor_locator(plt.NullLocator())
        ax1[pos[0]][pos[1]].xaxis.set_major_locator(plt.FixedLocator([1.0e-3,1.0e-5] ))
        ax1[pos[0]][pos[1]].xaxis.set_minor_locator(plt.NullLocator())
        ax1[pos[0]][pos[1]].margins(x=0)
####################################################################################
    if limit_plot:
        fig1.set_figheight(4.5)
        fig1.set_figwidth(8)
        fig1.supylabel("$F_2$",rotation="horizontal",x=0.02,y=0.95,fontsize=13)
        fig1.supxlabel("$x$",rotation="horizontal",y=0.015,x=0.99,fontsize=13)
        if saveflag:
            fig1.savefig(save1)
        else:
        	plt.show()
        return(0)
######################################################################################

    for i in range(q2len-(q2len//2+extra)):
        pos=[i//col2,(i)%col2]
        q2val=q2set[i+(q2len//2+extra)];
        #print(q2val);
        arr=[ ]
        arr2=[ ]
        #yarr=[]
        for j in range(len(q2)):
            if(q2[j]==q2set[i+(q2len//2+extra)]):
                arr.append([ x[j],f2[j] ])
                arr2.append([ x2[j],f22[j] ])
            
        arr=np.transpose(sorted(arr,key=sortingkey))
        arr2=np.transpose(sorted(arr2,key=sortingkey))
        #print( arr)
        xarr=[float(k) for k in arr[0]]
        yarr=[float(k) for k in arr[1]]
        xarr2=[float(k) for k in arr2[0]]
        yarr2=[float(k) for k in arr2[1]]
        #print( xarr, yarr)
        ax2[pos[0]][pos[1]].plot(xarr,yarr ,c='blue',ls="--")
        ax2[pos[0]][pos[1]].plot(xarr2,yarr2 ,c='red')
        #chi=0
        #count=0
        chi1.append([])
        chi2.append([])
        for j in range(len(q2d)):
            if(q2d[j]==q2set[i+q2len//2+extra]):
                ax2[pos[0]][pos[1]].errorbar(float(xd[j]),float(f2d[j]),yerr=float(f2e[j]),c='green')
                ax2[pos[0]][pos[1]].scatter(float(xd[j]),float(f2d[j]),s=12,marker="x",c='green')
                chi1[i+(q2len//2+extra)].append(pow( (float(f2d[j])-float(f2c[j]))/float(f2e[j]) ,2))
            if(q2d2[j]==q2set[i+q2len//2+extra]):
                chi2[i+(q2len//2+extra)].append(pow( (float(f2d2[j])-float(f2c2[j]))/float(f2e2[j]) ,2))
                
        ax2[pos[0]][pos[1]].set(xscale="log" ,   yscale='linear' )
        #ax2[pos[0]][pos[1]].legend()
        ax2[pos[0]][pos[1]].text(1.0e-4, 0.6,"$Q^2={val}\\;\\mathrm{{GeV^2}}$".format(val=float(q2val)),fontsize=7) 
        #ax2[pos[0]][pos[1]].yaxis.set_major_locator(plt.FixedLocator([0.01,0.05,0.1,0.5,1,2.5]))
        ax2[pos[0]][pos[1]].yaxis.set_major_locator(plt.LogLocator(base=10))
        ax2[pos[0]][pos[1]].xaxis.set_major_locator(plt.FixedLocator([1.0e-6,1.0e-4,1.0e-2] ))
        ax2[pos[0]][pos[1]].xaxis.set_minor_locator(plt.NullLocator())
        ax2[pos[0]][pos[1]].margins(x=0)
        
        
    #ax2[0][0].set(label="$F_2$",loc="top")
    #ax1[0][0].set(label="$F_2$",loc="top")
        
        #ax.scatter(x,f2,marker=".");
        #ax.set(  xscale="log" ,   yscale='linear' )
        #ax2[pos[0]][pos[1]].set_ylabel( "F_2" ,rotation="horizontal",loc='top')
    
    #ax1.set_ylabel( "F_2" ,rotation="horizontal",loc='top')
    fig1.set_figheight(8)
    fig1.set_figwidth(6)
    #fig.set_ylabel( "F_2" ,rotation="horizontal",loc='top')
    fig1.supylabel("$F_2$",rotation="horizontal",x=0.02,y=0.95,fontsize=13)
    fig2.supylabel("$F_2$",rotation="horizontal",x=0.02,y=0.95,fontsize=13)
    
    fig1.supxlabel("$x$",rotation="horizontal",y=0.015,x=0.99,fontsize=13)
    fig2.supxlabel("$x$",rotation="horizontal",y=0.015,x=0.99,fontsize=13)
    #fig1.subplots_adjust(bottom=0.05, right=0.95, top=0.95, left=0.05)
    #fig2.subplots_adjust(bottom=0.05, right=0.95, top=0.95, left=0.05)
    
    
    fig2.set_figheight(8)
    fig2.set_figwidth(6)
    #plt.yticks(y,[0.1,1])
    #fig1.tight_layout()
    #fig2.tight_layout()
    if saveflag:
        fig1.savefig(save1)
        fig2.savefig(save2)
    else:
    	plt.show()
        
    
    ################################## chi #########################################
    chipp=[sum(i) for i in (chi1)]
    print(sum([ sum(i) for i in chi1])/sum([ len(i) for i in chi1]))
    chipp2=[sum(i) for i in (chi2)]
    print(sum([ sum(i) for i in chi2])/sum([ len(i) for i in chi2]))
    
    diff=[chipp2[i]-chipp[i] for i in range(len(chipp)) ]
    
    chipp=[sum(i)/len(i) for i in (chi1)]
    chipp2=[sum(i)/len(i) for i in (chi2)]
    
    #chipp2=[sum(i)/len(i) for i in chi2]
    binlabel=[]
    for i in range(len(chipp)):
         if (i%(len(chipp)//4) ==0 ):
             binlabel.append("{val:.2e}".format(val=float(q2set[i]) ) )
         else:
             binlabel.append("")
    
    fig, ax1=plt.subplots(1,1,constrained_layout=True, sharex=True,sharey=True)
    #fig=plt.figure(constrained_layout=True)
    #ax3.clf()
    #ax1=fig.add_subplot(1,1,1)
    #ax2=fig.add_subplot(1,1,2,sharex=ax1,sharey=ax1)
    #ax3=fig.add_subplot(1,1,3,sharex=ax1)
    leg=[]
    leg.append(ax1.bar(q2set,chipp,fill=False,tick_label=binlabel,edgecolor="blue"))
    leg.append(ax1.bar(q2set,chipp2,fill=False,tick_label=binlabel,edgecolor="red"))
    #ax2.bar(q2set,chipp2,tick_label=binlabel,fill=False)
    
    ax1.grid(visible='true', axis='y')
    ax1.margins(x=0)
    ax1.legend(leg,["Without Sudakov", "With Sudakov"])
    #ax2.grid(visible='true', axis='y')
    #ax2.margins(x=0)
    
    
    
    #ax2.set_ylabel("$\chi^2/N$ ", loc='center' )
    #ax1.set_ylabel("$\chi^2/N$", loc='center' )
    ax1.set_ylabel("$\\chi^2/N$", loc='center' )
    #ax2.set_ylabel("With Sudakov\n$\\chi^2/N$", loc='center' )
    
    
    
    #fig.set_figheight(4)
    #fig.set_figwidth(12)
    
    #if saveflag:
    #    fig.savefig(save3)
    #else:
    #	plt.show()
    #############################################################
    #fig, ax3=plt.subplots(1,1,sharey=True,sharex=True,constrained_layout=True)
    #ax3.bar(q2set,diff,tick_label=binlabel,color='purple',edgecolor="purple",fill=True)
    #ax3.set_xlabel("$Q^2\\;(\\mathrm{GeV^2})$", loc='right' )
    #ax3.margins(x=0)
    #ax3.grid(visible='true', axis='y')
    #ax3.set_ylabel("difference\n$\\chi^2$", loc='center' )
    fig.set_figheight(2)
    fig.set_figwidth(8)
    
    if saveflag:
        fig.savefig(save3)
    else:
    	plt.show()
    
    return(0)
main()




















