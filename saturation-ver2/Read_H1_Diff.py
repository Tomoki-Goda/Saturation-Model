#! /usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys, getopt , os

def main():

    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:],"i:o:",["in","out"] )
    except getopt.GetoptError as err:
        print("Error reading option")

    for opt, arg in opts:
        if opt in ['-i','-in']:
            read_file=arg
        elif opt in ['-o','-out']:
            outdir=arg


    #read_file="./data/H1diff2006.txt"
    data=pd.read_csv(read_file,sep=" ",index_col=None)
    print(data)
    #print(sorted(list(set(data["Q2"].values))))
    data= data[data['Q2']<=90]
    if not('xp' in data):
        xp=data['x']/data['beta']
        data['xp']=xp

    Q2_set=sorted(list(set(data["Q2"].values)))
    print(sorted(list(set(data["beta"].values))))
    print(sorted(list(set(data["xp"].values))))
    #input()
    print("Q2 ",Q2_set,"\n\n")
    for Q2 in Q2_set:
        data_Q2=data[data['Q2']==Q2]
        beta_set=sorted(list(set(data_Q2["beta"].values)))
        print("beta  ",  beta_set,"\n\n")
        for beta in beta_set:
            data_beta=data_Q2[data_Q2['beta']==beta]
            #xp_set=sorted(list(set(data_beta["xp"].values)))
            xp_set=sorted(list(set(data_beta["xp"].values)))
            print(data_beta)
            xmax=max(xp_set)
            xmin=min(xp_set)

            os.system("{DIR}/diffraction -in {DIR}/result.txt -beta {beta_:} -Q2 {q2_:} -xmin {xmin_:} -xmax {xmax_:} -out {DIR}/file-{q2_:}-{beta_:}.txt ".format(DIR=outdir,q2_=Q2,beta_=beta,xmax_=xmax,xmin_=xmin ))
            fd3=[]
            xp=[]
            with open("{DIR}/file-{q2_:}-{beta_:}.txt".format(q2_=Q2,beta_=beta,DIR=outdir)) as fi:
                for line in fi:
                    plot_data=line.strip().split("\t")
                    fd3.append(float(plot_data[1]))
                    xp.append(float(plot_data[0]))

            plt.plot(xp,fd3)
            #datapt=np.transpose(frame[['xp','xF']].values)
            datapt=np.transpose(data_beta[['xp','xp*sigmaD(3)']].values)
            error=data_beta['dtot'].values
            error=error*datapt[1]/100
            #error=frame[["stat","sys+","sys-"]].values
            #errp=[np.sqrt(float(i[0])**2+float(i[1])**2 ) for i in error]
            #errn=[np.sqrt(float(i[0])**2+float(i[2])**2) for i in error]
            #plt.errorbar([float(i) for i in datapt[0]],[float(i) for i in datapt[1]],yerr=[errn,errp] , c='b')
            plt.errorbar([float(i) for i in datapt[0]],[float(i) for i in datapt[1]],yerr=error , c='b')
            plt.scatter([float(i) for i in datapt[0]],[float(i) for i in datapt[1]] ,c='black')
            
            plt.xscale('log')
            
            #plt.show()
            plt.savefig("{DIR}/plot-{q2_:}-{beta_:}.png".format(DIR=outdir,beta_=beta,q2_=Q2))
            plt.clf()


    

if __name__=="__main__":
    main()




