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
    label=["" for i in range(len(args))] 
    for opt,arg in opts:
        if opt in ["-l","-label"]:
            label=arg.split(":")
           

    data=[]
    with open(args[0],"r") as fi:
        for line in fi:
            row=line.strip().split("\t")
            data.append(row)
    data=pd.DataFrame(data,columns=["x","Q2","data","err","val"])
    comp=[]
    
    for i in range(1,len(args)):
        data2=[]
        with open(args[i],"r") as fi:
            for line in fi:
                row=line.strip().split("\t")
                data2.append(row)
        comp.append(pd.DataFrame(data2,columns=["x","Q2","val"]))
    
    Q2set=data['Q2'].values
    Q2set=set(Q2set)
    Q2set=sorted(list(Q2set),key=float)
    length=len(Q2set)
    col=4
    row=length//(2*col)
    if row*2*col<length:
        row+=1
    fig1, ax1= plt.subplots(row,col,sharey='row',sharex='col',layout="constrained")
    fig2, ax2= plt.subplots(row,col,sharey=True,sharex=True,layout="constrained")
   
    for i,Q2 in enumerate(Q2set): 
        frame=data.where(data['Q2']==Q2).dropna().astype(float)
        frame2=[i.where(i['Q2']==Q2).dropna().astype(float) for i in comp]
        
        frame=frame.sort_values("x")
        frame2=[i.sort_values("x") for i in frame2]
        x=[float(i) for i in frame["x"].tolist()]
        y=[float(i) for i in frame["data"].tolist()]
        e=[float(i) for i in frame["err"].tolist()]
        x2=[[float(i) for i in j["x"].tolist()] for j in frame2]
        y2=[[float(i) for i in j["val"].tolist()] for j in frame2]
        
            
        if i//4<row:
            #print(1,": ", i," ",i//4,"  ",i%4)
            ax1[i//col][i%col].errorbar(x,y,yerr=e ,c='green',ls='none' )
            ax1[i//col][i%col].scatter(x,y,s=5,marker="x",c='green' )
            for k in range(len(x2)):
                ax1[i//col][i%col].plot(x2[k],y2[k],c='b',linewidth=1 )
                
            ax1[i//col][i%col].set(xscale='log')
            ax1[i//col][i%col].text(0.05,0.9,"$Q^2={val}\\;\\mathrm{{GeV^2}}$".format(val=float(Q2)),fontsize=9,transform=ax1[i//col][i%col].transAxes)
            ax1[i//col][i%col].xaxis.set_major_locator(plt.FixedLocator([1.0e-3,1.0e-5] ))
        else:
            if i>length-1:
                break
            #print(2,": ", i,"  ",(i-4*row)//4,"  ",(i-4*row)%4)
            ax2[(i-col*row)//col][(i-col*row)%col].errorbar(x,y,yerr=e ,c='green' ,ls='none' )
            ax2[(i-col*row)//col][(i-col*row)%col].scatter(x,y, s=5,marker="x",c='green')
            leg=[]
            for k in range(len(x2)):
                 leg.append(ax2[(i-col*row)//col][(i-col*row)%col].plot(x2[k],y2[k],c='b',linewidth=1 ))
            ax2[(i-col*row)//col][(i-col*row)%col].set(xscale='log')
            ax2[(i-col*row)//col][(i-col*row)%col].text(0.05,0.9,"$Q^2={val}\\;\\mathrm{{GeV^2}}$".format(val=float(Q2)),fontsize=9,transform=ax2[(i-col*row)//col][(i-col*row)%col].transAxes) 
            ax2[(i-col*row)//col][(i-col*row)%col].xaxis.set_major_locator(plt.FixedLocator([1.0e-3,1.0e-5] ))
        
    
    ax2[row-1][col-1].legend(leg[:][0],label)
    fig1.get_layout_engine().set(hspace=None, wspace=None)
    fig2.get_layout_engine().set(hspace=None, wspace=None)
    fig1.set_figheight(8)
    fig1.set_figwidth(6)
    fig1.supylabel("$F_2$",rotation="horizontal",x=0.02,y=0.98,fontsize=13)
    fig2.supylabel("$F_2$",rotation="horizontal",x=0.02,y=0.98,fontsize=13)
    
    fig1.supxlabel("$x$",rotation="horizontal",y=0.015,x=0.99,fontsize=13)
    fig2.supxlabel("$x$",rotation="horizontal",y=0.015,x=0.99,fontsize=13)
    
    
    fig2.set_figheight(8)
    fig2.set_figwidth(6)
    plt.show()
        
        
        
        
    
    
if __name__=="__main__":
    main()
