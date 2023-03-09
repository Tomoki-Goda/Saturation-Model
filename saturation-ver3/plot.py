#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

import math
import sys 
import getopt

#dataarray=np.array([]);
def import_array(name ):
	dataarray=[];
	with open(name,"r") as fi:
	#with open("./pwf.txt","r") as fi:
		for line in fi:
			#data=line.readline();
			data=line.strip().split('\t')
			#data=[math.log(float(j))/math.log(10) for j in data]
			data=[float(j) for j in data]
			dataarray.append(data)
	dataarray=np.array(dataarray)
	dataarray=np.transpose(dataarray)
	return(dataarray)

def main():
    xs="log"
    ys="log"
    sflag=False
    sname=""
    axes=["",""]
    plottitle=""
    column=[0,1]
    xlim=[0,0]
    ylim=[0,0]
    argv=sys.argv[1:]
    plotcolour=["red", "blue", "green", "cyan", "magenta", "brown", "orange", "purple", "yellow"]
    scatter=False
    labels=[]
 
    fig,ax=plt.subplots(1,1,layout="constrained")
 
    try:
        opts , args =getopt.getopt(argv,"x:y:l:s:a:t:p:c:n:X:Y:S" )
        labels=args;
        plotstyle=["-" for i in range(len(args))]
        

    except:
        print("option eror");
    
    for opt,arg in opts:
        if opt in ["-x","--xscale"]:
            xs=arg;
        elif opt in ["-y", "--yscale"]:
            ys=arg;
        elif opt in ["-l","--labels"]:
            labels=arg.split(":")
            print(labels)
        elif opt in ["-s", "--save"]:
            sflag=True
            sname=arg
        elif opt in ["-a","--axes"]:
            axes=arg.split(":");
        elif opt in ["-t" , "--title"]:
            plottitle=arg
        elif opt in ["-p","--plot-style"]:
            plotstyle=arg.split();
        elif opt in ["-c", "--colour"]:
            plotcolour=arg.split()
        elif opt in ["-X","-Xlim"]:
        	xlim=arg.split()
        	xlim=[float(i) for i in xlim]
        elif opt in ["-Y","-Ylim"]:
        	ylim=arg.split()
        	ylim=[float(i) for i in ylim]
        elif opt in ["-S","--scatter"]:
            scatter=True
        elif opt in ["-n","-number"]:
            column=[int(i) for i in arg.split()];
            


    data_array=[];
    counter=0;
    legs=[];
    for i in args:
        print(i)
        data_array=import_array(i)
        if scatter:
            leg=ax.scatter(data_array[column[0]],data_array[column[1]],label=labels[counter],s=1 ,color=plotcolour[counter] ) 
        else:
            leg=ax.plot(data_array[column[0]],data_array[column[1]],label=labels[counter],linestyle=plotstyle[counter] ,color=plotcolour[counter] ) 
        #if counter<len(labels):
        if (labels!=[] and labels[counter]!=""):
        	legs.append(leg[0])
        counter+=1
    print("LABELS ", labels);
    if "" in labels:
        labels.remove("")
        if len(legs)>0:
            ax.legend(legs,labels)
    ax.set(title=plottitle,  ylabel=axes[1],    xlabel=axes[0],  xscale= xs ,   yscale=ys )
    if(xlim[0]!=xlim[1]):
        ax.set_xlim(xlim[0],xlim[1])
    
    if(ylim[0]!=ylim[1]):
        ax.set_ylim(ylim[0],ylim[1])
    	
    ax.grid("true")
    ax.legend()
   
    
    
    fig.set_figheight(4);
    fig.set_figwidth(5);

    if sflag:
        fig.savefig(sname)
    else:
        plt.show()

main()




