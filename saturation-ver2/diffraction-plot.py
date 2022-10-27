#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import getopt , sys
import pandas as pd
import os

def main():
    try:
        opts, args=getopt.gnu_getopt(sys.argv[1:],"i:o:",["in","out"])
    except getopt.GetoptError as err:
        print("error")

    print(opts);
    for opt,arg in opts:
        if opt in ["-i","-in"]:
            input_file=arg
            print("in ", arg);
        elif opt in ["-o","-out"]:
            output_file=arg

    data=[]
    with open(input_file,"r") as fi:
        for line in fi:
            row=line.strip().split("\t")
            data.append(row)
    data=pd.DataFrame(data,columns=["Q2","beta",'xp','xF', 'stat','sys+','sys-' ])
    #print(data)
    Q2set=data['Q2'].values
    Q2set=set(Q2set)
    Q2set=np.array(list(Q2set))
    Q2set.sort()
    for q2 in Q2set:
        tf=data["Q2"]==q2
        betaframe=data.where(tf).dropna()
        betaset=np.array(list(set(betaframe['beta'].values)))
        betaset.sort()
        #print(betaset)
        for beta in betaset:
            tf2=data["beta"]==beta
            frame=betaframe.where(tf2).dropna()
            print(frame)
            betaset=np.array(list(set(frame['xp'].values)))
            xmax=max(betaset)
            xmin=min(betaset)

            os.system("echo file -in reufile -beta {beta_:} -Q2 {q2_:} -xmin {xmin_:} -xmax {xmax_:} -out file-{q2_:}-{beta_:}.txt ".format(q2_=q2,beta_=beta,xmax_=xmax,xmin_=xmin ))

if __name__=="__main__":
    main()
















