#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import getopt , sys
import pandas as pd
import os

def main():
    try:
        opts, args=getopt.gnu_getopt(sys.argv[1:],"",[])
    except getopt.GetoptError as err:
        print("error")

    print(opts);
    #for opt,arg in opts:
        #if opt in ["-i","-in"]:
        #    input_file=arg
        #   print("in ", arg);
        #elif opt in ["-o","-out"]:
        #   outdir=arg

    data=[]
    with open(args[0],"r") as fi:
        for line in fi:
            row=line.strip().split("\t")
            data.append(row)
    data=pd.DataFrame(data,columns=["x","Q2","data","err","val"])
    print(data)
    Q2set=data['Q2'].values
    Q2set=set(Q2set)
    Q2set=sorted(list(Q2set),key=float)
    
    length=len(Q2set)
    row=length//(2*4)
    if row*8<length:
    	row+=1
    fig1, ax1= plt.subplots(row,4,sharey=True,sharex=True,layout="constrained")
    fig2, ax2= plt.subplots(row,4,sharey=True,sharex=True,layout="constrained")
    for i,Q2 in enumerate(Q2set): 
        frame=data.where(data['Q2']==Q2).dropna().astype(float)
        
        frame=frame.sort_values("x")
        #print(frame)
        #x1=[[float(j) for j in i] for i in frame["x",'data'].tolist()]
        
        x=[float(i) for i in frame["x"].tolist()]
        y=[float(i) for i in frame["data"].tolist()]
        e=[float(i) for i in frame["err"].tolist()]
        
        x2=[float(i) for i in frame["x"].tolist()]
        y2=[float(i) for i in frame["val"].tolist()]
        
        #for k,yk in enumerate(y):
        leng=len(x)
        for k in range(leng) :
           if y[leng-k-1]==0.0:
                x.pop(leng-k-1)
                y.pop(leng-k-1)
                e.pop(leng-k-1)
        
            
        if i//4<row:
            #print(1,": ", i," ",i//4,"  ",i%4)
            ax1[i//4][i%4].errorbar(x,y,yerr=e ,c='green',ls='none' )
            ax1[i//4][i%4].scatter(x,y,s=5,marker="x",c='green' )
            ax1[i//4][i%4].plot(x2,y2,c='b' )
            ax1[i//4][i%4].set(xscale='log')
            ax1[i//4][i%4].text(0.05,0.9,"$Q^2={val}\\;\\mathrm{{GeV^2}}$".format(val=float(Q2)),fontsize=9,transform=ax1[i//4][i%4].transAxes) 
        else:
            if i>length-1:
            	break
            #print(2,": ", i,"  ",(i-4*row)//4,"  ",(i-4*row)%4)
            ax2[(i-4*row)//4][(i-4*row)%4].errorbar(x,y,yerr=e ,c='green' ,ls='none' )
            ax2[(i-4*row)//4][(i-4*row)%4].scatter(x,y, s=5,marker="x",c='green')
            ax2[(i-4*row)//4][(i-4*row)%4].plot(x2,y2,c='b' )
            ax2[(i-4*row)//4][(i-4*row)%4].set(xscale='log')
            ax2[(i-4*row)//4][(i-4*row)%4].text(0.05,0.9,"$Q^2={val}\\;\\mathrm{{GeV^2}}$".format(val=float(Q2)),fontsize=9,transform=ax2[(i-4*row)//4][(i-4*row)%4].transAxes) 
        
    fig1.get_layout_engine().set(hspace=None, wspace=None)
    fig2.get_layout_engine().set(hspace=None, wspace=None)
    fig1.set_figheight(8)
    fig1.set_figwidth(6)
    fig1.supylabel("$F_2$",rotation="horizontal",x=0.02,y=0.95,fontsize=13)
    fig2.supylabel("$F_2$",rotation="horizontal",x=0.02,y=0.95,fontsize=13)
    
    fig1.supxlabel("$x$",rotation="horizontal",y=0.015,x=0.99,fontsize=13)
    fig2.supxlabel("$x$",rotation="horizontal",y=0.015,x=0.99,fontsize=13)
    
    
    fig2.set_figheight(8)
    fig2.set_figwidth(6)
    plt.show()
        
        
        
        
    
    
if __name__=="__main__":
    main()
