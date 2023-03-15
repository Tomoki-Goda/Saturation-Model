#! /usr/bin/env python3

import numpy as np
import getopt,sys

def main():
    
    try:
       opts, args = getopt.gnu_getopt(sys.argv[1:],"")
              

    except:
        print("option eror");
    
    #for opt,arg in opts:
       
            
    name=args[0]
    name2=args[1]
    name3=args[2]
    daraarray=[]
    dataarray2=[]
    
    with open(name,"r") as fi:
        with open(name2,"r") as fi2:
            with open(name3,"w") as fi3:
                for line in fi:
                    line2=fi2.readline()
                    data=line.strip().split('\t')
                    data2=line2.strip().split('\t')
                    #print("{0}\t{1}\t{2}".format(data[0],data[1],float(data[2])-float(data2[2])))
                    assert data[0]==data2[0]
                    assert data[1]==data2[1]
                    fi3.write("{0}\t{1}\t{2}\n".format(data[0],data[1],float(data[2])-float(data2[2])))
                    #print("{0}\t{1}\t{2}".format(data[0],data[1],float(data[2])-float(data2[2])))
if __name__=="__main__":
    main()
