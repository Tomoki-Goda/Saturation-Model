#! /usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys, getopt , os

def main():
    plot=True
    compute=False
    parallelkernel=4;
    indir=["."]
    savedir="."
    outdir="."
    readfile=""
    plotcolor=['b','r','g','m']
    plotstyle=['--','-',':','-.']
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:],"d:hi:o:ps:cj:",['data',"help","in","out","plot",'save','compute','j'] )
    except getopt.GetoptError as err:
        print("Error reading option")

    for opt, arg in opts:
        if opt in ['-i','-in']:
            rundir=arg
        elif opt in ['-d','-data']:
            read_file=arg
        elif opt in ['-o','-out']:
            outdir=arg
            if compute:
                indir=[outdir]
        elif opt in ['-p','-plot']:
            plot=True
            compute=False
            #indir=arg.strip().split(':');
        elif opt in ['-s','--save']:
            savedir=arg
        elif opt in ['-c','--compute']:
            compute=True
        elif opt in ['-j','--j']:
            parallelkernel=arg
        elif opt in ['-h','--help']:
            print(' getopt.gnu_getopt(sys.argv[1:],"h:i:o:ps:c",["help","in","out","plot","save","compute"]' )
            return(0)

    if plot:
        indir=args

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
    lengthQ2=len(list(set(data['Q2'])))
    lengthbeta=len(list(set(data['beta'])))

    #plt.ion()
    if plot:
        plt.margins(0)
        fig,ax=plt.subplots(lengthQ2,lengthbeta, sharex=True,sharey=True,layout="constrained")
        fig.get_layout_engine().set(wspace=0,hspace=0,w_pad=0,h_pad=0)
        #fig.xscale('log')
        print(lengthQ2, ",  " ,lengthbeta)
        fig.set_figheight(lengthQ2)
        fig.set_figwidth(1.2*lengthbeta)
    beta_set=sorted(list(set(data["beta"].values)))
    
    parallelargs=[]
    parallelargs_All=[]
    args=[]

    print("Q2 ",Q2_set,"\n\n")
    for i in range(len(Q2_set)):
        Q2=Q2_set[i]
        data_Q2=data[data['Q2']==Q2]
        #beta_set=sorted(list(set(data_Q2["beta"].values)))
        print("beta ",i," ",  beta_set,"\n\n")
        #parallelargs=[]
        for j in range(len(beta_set)):
            beta=beta_set[j]
            data_beta=data_Q2[data_Q2['beta']==beta]
            if len(data_beta)==0:
                continue
            #xp_set=sorted(list(set(data_beta["xp"].values)))
            x_set=sorted(list(set(data_beta["x"].values)))
            #print(data_beta)
            xmax=max(x_set)
            xmin=min(x_set)
            if compute:
                #os.system("{DIR}/diffraction -in {DIR}/result.txt -beta {beta_:} -Q2 {q2_:} -xmin {xmin_:} -xmax {xmax_:} -out {DIR}/file-{q2_:}-{beta_:}.txt ".format(DIR=outdir,q2_=Q2,beta_=beta,xmax_=xmax,xmin_=xmin ))
                parallelargs.append([rundir+"/result.txt","{dir_}/file-{q2_:}-{beta_:}.txt".format(dir_=outdir,q2_=Q2,beta_=beta),'{}'.format(Q2),'{}'.format(beta),'{}'.format(xmin),'{}'.format(xmax)])

            else:
                for k in range(len(indir)):
                    fd3=[]
                    xp=[]
                    if(os.path.exists("{DIR}/file-{q2_:}-{beta_:}.txt".format(q2_=Q2,beta_=beta,DIR=indir[k])) ):
                        with open("{DIR}/file-{q2_:}-{beta_:}.txt".format(q2_=Q2,beta_=beta,DIR=indir[k])) as fi:
                            for line in fi:
                                plot_data=line.strip().split("\t")
                                fd3.append(1.23*float(plot_data[1]))
                                xp.append(float(plot_data[0]))
                        ax[i][j].plot(xp,fd3,c=plotcolor[k],ls=plotstyle[k])
                        #datapt=np.transpose(frame[['xp','xF']].values)
                        datapt=np.transpose(data_beta[['xp','xp*sigmaD(3)']].values)
                        error=data_beta['dtot'].values
                        error=error*datapt[1]/100
                        #error=frame[["stat","sys+","sys-"]].values
                        #errp=[np.sqrt(float(i[0])**2+float(i[1])**2 ) for i in error]
                        #errn=[np.sqrt(float(i[0])**2+float(i[2])**2) for i in error]
                        #plt.errorbar([float(i) for i in datapt[0]],[float(i) for i in datapt[1]],yerr=[errn,errp] , c='b')
                    ax[i][j].errorbar([float(i) for i in datapt[0]],[float(i) for i in datapt[1]],yerr=error , c='black',marker="x",fmt='none')
                ax[i][j].scatter([float(i) for i in datapt[0]],[float(i) for i in datapt[1]] ,c='black',s=2)
                    
                ax[i][j].set(xscale='log')
                ax[i][j].tick_params(direction='in')
                #plt.show()i
        #if compute:
        #    parallelargs_All.append(parallelargs)
    
    if compute:
        #paralellargs_All=np.transpose(parallelargs_All);
        if parallelkernel=="1":
            for i in parallelargs:
                command=rundir+"/diffraction "
                for j in i:
                    command=command+j+" "
                print(command);
                os.system(command)

        else: 
            parallelargs=np.transpose(parallelargs)
            command="parallel  -j "+str(parallelkernel)+" --link "
            command+=rundir+"/diffraction "

            for i in parallelargs:
                command=command+"::: "
                for j in i:
                    command=command+j+" "
            #print(command)
            os.system(command)

    if plot:
        plt.savefig("{DIR}/plot-diff.png".format(DIR=savedir,beta_=beta,q2_=Q2))
    #plt.clf()


    

if __name__=="__main__":
    main()




