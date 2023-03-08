#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import getopt , sys
import pandas as pd
import os

def main():
    try:
        opts, args=getopt.gnu_getopt(sys.argv[1:],"l:p:c:h",["label","plot-style","color",'help'])
    except getopt.GetoptError as err:
        print("error")

    print(opts);
    label=["" for i in range(len(args))] 
    color=['b','r','g']
    ls=['-.','-.','-.']
    for opt,arg in opts:
        if opt in ["-l","-label"]:
            label=arg.split(":")
        elif opt in ["-c","-color"]:
            color=arg.split(":")
        elif opt in ["-p","-plot-style"]:
            ls=arg.split(" ")
        elif opt in ['-h','--help']:
            print("-l label (: separated)\n-c color (: separated)\n-p plot-style (space separated)\n args=F2_exp_data...\n ")
            exit(0)
           

    
 
    comp=[]
    
    for i in range(len(args)):
        data2=[]
        with open(args[i],"r") as fi:
            for line in fi:
                row=line.strip().split("\t")
                data2.append(row)
        comp.append(pd.DataFrame(data2,columns=["x","Q2","data","err","val"]))
    data=comp[0]
    Q2set=data['Q2'].values
    Q2set=set(Q2set)
    Q2set=sorted(list(Q2set),key=float)
    length=len(Q2set)
    col=4
    row=length//(2*col)
    if row*2*col<length:
        row+=1
    fig1, ax1= plt.subplots(1,1,layout="constrained")
    binlabel=[]
    for i,j in enumerate(Q2set):
         if (i%(len(Q2set)//4) ==0 ):
             binlabel.append("{val:.2e}".format(val=float(j) ) )
         else:
             binlabel.append("")
    #fig1.set_constrained_layout_pads(hspace=0,wspace=0)
    hisall=[]
    for i,Q2 in enumerate(Q2set): 
        #frame=data.where(data['Q2']==Q2).dropna().astype(float)
        frame=[i.where(i['Q2']==Q2).dropna().astype(float) for i in comp]
        
        his=[]
        for j in frame:
            chisum=0
            d1=[float(i) for i in j["data"].tolist()]
            d2=[float(i) for i in j["err"].tolist()]
            d3=[float(i) for i in j["val"].tolist()]
            for k in range(len(d1)):
                chisum+=pow((d1[k]-d3[k])/d2[k],2)
            print(chisum)
            his.append(chisum/len(d1))
        hisall.append(his)
    hisall=np.transpose(hisall)
    leg=[]
    for i,j in enumerate(hisall):
    	leg.append(ax1.bar(Q2set,j,fill=False,tick_label=binlabel,edgecolor=color[i],ls=ls[i]))
    fig1.get_layout_engine().set(hspace=None, wspace=None)
    fig1.set_figheight(2)
    fig1.set_figwidth(8)
    
    ax1.grid(visible='true', axis='y')
    ax1.margins(x=0)
    ax1.legend(np.transpose(leg)[0],label)
    
    ax1.set_ylabel("$\\chi^2/N$", loc='center' )
    plt.show()
        
        
        
        
    
    
if __name__=="__main__":
    main()
