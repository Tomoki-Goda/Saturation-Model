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
        
    with open(args[0]+"/data.txt", "r") as fi:
    	data1=[]
    	for i in fi:
    		data1.append(i.strip().split('\t'))
     with open(args[1]+"/data.txt", "r") as fi:
    	data2=[]
    	for i in fi:
    		data2.append(i.strip().split('\t'))
    	
    		
