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

    argv=sys.argv[1:]
    plotcolour=["red", "blue", "green", "cyan", "magenta", "brown", "orange", "purple", "yellow"]
 	
 
    fig,ax=plt.subplots()
 
    try:
        opts , args =getopt.getopt(argv,"x:y:l:s:a:t:p:c:" )
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

    data_array=[];
    counter=0;
    for i in args:
        print(i)
        data_array=import_array(i)
        ax.plot(data_array[0],data_array[1],label=labels[counter],linestyle=plotstyle[counter] ,color=plotcolour[counter] ) 
        counter+=1
    
    ax.set(title=plottitle,  ylabel=axes[1],    xlabel=axes[0],  xscale= xs ,   yscale=ys )
    ax.grid("true")
    ax.legend()
    if sflag:
        fig.savefig(sname)

    plt.show()

main()




