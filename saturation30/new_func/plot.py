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

    argv=sys.argv[1:]
    

    try:
        opts , args =getopt.getopt(argv,"x:y:" )

    except:
        print("option eror");
    
    for opt,arg in opts:
        if opt in ["-x","--xscale"]:
            xs=arg;
        elif opt in ["-y", "--yscale"]:
            ys=arg;

    data_array=[];
    for i in args:
        print(i)
        data_array=import_array(i)
        plt.plot(data_array[0],data_array[1],label=i )
    
    plt.xscale(xs)
    plt.yscale(ys)
    plt.legend()
    plt.show()

main()




