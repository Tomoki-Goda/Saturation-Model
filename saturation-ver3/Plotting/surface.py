#! /usr/bin/env python3
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
import getopt,sys



#fig, ax = plt.subplots()

def main():
    grids="N/A"
    posQ2=0
    savefile=""
    dipole=False
    contour=False
    run=False
    alpha=1
    rapidity=False
    logz=False
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "Yrdmhs:g:q:cl", ['Y','run','dipole','multiple','help',"save=","grid=","Q2=","contour","log"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
    for opt, arg in opts:
        print(opt," ")
        if opt in ["-s","--save"]:
            savefile=arg
        elif opt in ["-g","--grid"]:
            grids=[arg]
        elif opt in ["-q","--Q2"]:
            print(arg)
            posQ2=int(arg)
        elif opt in ['-d','--dipole']:
            dipole=True
        elif opt in ['-h', '--hellp']:
            print("-s <save file>\n-g <grid>\n-q <position of Q2 in the list of Q2 in the grid...> ")
            print("-m to say multiple grid files, and list grids files after options.\nuse -d to indicate it is  dipole cross section. ")
            sys.exit()
        elif opt in ['-m','--multiple']:
            grids=args;
        elif opt in ['-c','--contour']:
            contour=True
        elif opt in ['-r','--run']:
            run=True
        elif opt in ['-Y','--Y']:
            rapidity=True
        elif opt in ['-l', '--log']:
            logz=True
    if grids=="N/A":
        grids=args[0]

    #DIR="../GKSgluon-1.1/grids/"
    #with open("./KSgluon-2.0/grids/KSnonlinear.dat","r") as fi:
    #with open("../Run/GBW/Mass0.0-Qup650-Model0-Sud0/gluon-grid.dat","r") as fi:
    #with open("../Run/GBWS-Fix-S/Mass0.0-Qup650-Model22-Sud1/gluon-grid.dat","r") as fi:
    #with open(DIR+"KSnonlinear.dat","r") as fi:
    if contour:
        fig, ax =plt.subplots(1,len(grids), sharex=True,sharey=True,layout='constrained')
    else:
        fig, ax =plt.subplots(1,len(grids),layout='constrained',subplot_kw={'projection':'3d'})
    
    
    fig.set_constrained_layout_pads(hspace=0,wspace=0,h_pad=0.025,w_pad=0.025)
    fig.set_figheight(5)
    fig.set_figwidth(6*(len(grids)+0.1))
    
    if len(grids)==1:
        ax=[ax]
    for gpos in range(len(grids)):
        grid=grids[gpos]
        with open(grid,"r") as fi:
            X=[]
            Y=[]
            z=[]
            Z=[]
            x="n"
            counter=0
            lineno=0
            prevQ2=0
            Q2=set([])
            printflag=True
            for i in fi:
                data=i.strip().split()
                if len(data)==0:
                    continue
                if len(data)==4:
                    Q2.add(data[2])
                    if float(prevQ2)>float(data[2]):
                        lineno=0;
                    lineno+=1;
                    prevQ2=data[2]
                    if(posQ2!=(lineno-1)):
                        continue
                    elif printflag:
                        printflag=False
                        print("$Q^2$=",np.exp(float(data[2])));
                kt2=float(data[1])
                if run:
                    alpha=4*3.1415/(9*np.log((kt2+1)/0.09));
                     
                if((x=="n") or (data[0]!=x)):
                    if rapidity:
                        X.append(np.log(1/pow(10,float(data[0]))))
                    else:
                        X.append(float(data[0]))
                    if x!="n":
                        Z.append(z)
                        z=[]
                        counter+=1
                    z.append(alpha* float(data[len(data)-1]) )
                    x=data[0];
                else:            
                    z.append(alpha* float(data[len(data)-1]) )
                if(counter==0):
                    Y.append(kt2)
            Z.append(z)
        if len(Q2)>0:
            Q2=[np.log10(np.exp(float(i))) for i in list(Q2)]
            Q2.sort()
            print(Q2[0],"  ", Q2[len(Q2)-1]);
            
        if logz:
            Zn=[]
            for i, zi in enumerate(Z):
                zn=[]
                for j, zij in enumerate(zi):
                    if zij>0:
                        Z[i][j]=np.log10(zij)
                        zn.append(-50)
                    else:
                        zn.append(np.log10(-zij))
                        Z[i][j]=-50
                Zn.append(zn)
        #X=[np.log10(x) for x in X] 
        #Y=[np.log10(y) for y in Y] 
        print(X[0],"  ", X[len(X)-1]);
        print(Y[0],"  ", Y[len(Y)-1]);
        
        X, Y = np.meshgrid(X,Y )
        
        #ax[gpos].set(yscale="log",xscale="log")
        if logz:
            offset=-10
        else:
            offset=-0.01
        if dipole :
            
            lev=[i*0.1 for i in range(0,12)]
        else:
            lev=[(i-0.25)*0.00325 for i in range(-1,11)]
        if logz:
            #print(lev)
            lev=[(-10 if i<=0.0 else np.log10(i) )  for i in lev]
            print(lev)            
            #lev=list(set(lev))
        if contour:
            
            surf = ax[gpos].contourf(np.array(X),np.array(Y),np.transpose(np.array(Z)),levels=25,cmap=cm.coolwarm)
        else:
            ax[gpos].plot_wireframe(np.array(X),np.array(Y),np.transpose(np.array(Z)), rcount=4, ccount=0,color="r",ls="-." )
            ax[gpos].plot_wireframe(np.array(X),np.array(Y),np.transpose(np.array(Z)), rcount=0, ccount=4,color="b",ls="-")
            if dipole:
                surf = ax[gpos].contourf(np.array(X),np.array(Y),np.transpose(np.array(Z)),zdir='z',offset=offset,levels=lev,cmap=cm.cool)#,alpha=0.5)
                #surf = ax[gpos].contour(np.array(X),np.array(Y),np.transpose(np.array(Z)),zdir='z',offset=offset,levels=lev,colors='black',linewidths=0.5)
                #surf = ax[gpos].contour(surf,levels=surf.levels[::2],colors='black',linewidths=0.5)
            else:
                surf = ax[gpos].contourf(np.array(X),np.array(Y),np.transpose(np.array(Z)),zdir='z',offset=offset,levels=lev,cmap=cm.cool)#,alpha=0.5)
                #surf = ax[gpos].contour(np.array(X),np.array(Y),np.transpose(np.array(Z)),zdir='z',offset=offset,levels=5,colors='black',linewidths=0.5)
                #surf = ax[gpos].contour(surf,levels=surf.levels[::2],colors='black',linewidths=0.5)
            #surf = ax[gpos].contour(np.array(X),np.array(Y),np.transpose(np.array(Z)),zdir='z',offset=-0.01,levels=10,,colors='black')
            if logz:
                ax[gpos].plot_wireframe(np.array(X),np.array(Y),np.transpose(np.array(Zn)),rcount=4, ccount=0,color="g",ls="-." )
                ax[gpos].plot_wireframe(np.array(X),np.array(Y),np.transpose(np.array(Zn)), rcount=0, ccount=4,color="y",ls="-")
            if dipole:
                if logz:
                     ax[gpos].set(zlim=(offset,0))
                     ax[gpos].set_zlabel("$\\log_{10}(\\sigma(x,r)/\\sigma_0)$",fontsize=15)
                else:
                     ax[gpos].set(zlim=(offset,1.01))
                     ax[gpos].set_zlabel("$\\sigma(x,r)/\\sigma_0$",fontsize=15)
                ax[gpos].view_init(28, -28)
                #ax[gpos].set_ylabel('$\ln (r^2\\;[\\mathrm{GeV^2}])$',rotation='vertical',loc='top')
            else:
                if logz:
                    ax[gpos].set(zlim=(-25,-1.5))
                    ax[gpos].set_zlabel("$\\log_{10}(\\mathcal{F}(x,k^2)/\\sigma_0)$",fontsize=15)
                else:
                    ax[gpos].set(zlim=(-0.01,0.03))
                    ax[gpos].set_zlabel("$\\mathcal{F}(x,k^2)/\\sigma_0$",fontsize=15)
                ax[gpos].view_init(30, 25)
                #ax[gpos].set_ylabel('$\ln (k^2\\;[\\mathrm{GeV^2}])$',rotation='vertical',loc='top')
            ax[gpos].set_xlabel('$\\ln x$',loc='right',fontsize=15)
            if dipole:
                ax[gpos].set_ylabel('$\\ln (r\\;[\\mathrm{GeV^2}])$',rotation='vertical',loc='top',fontsize=15)
            else:
                ax[gpos].set_ylabel('$\\ln (k^2\\;[\\mathrm{GeV^2}])$',rotation='vertical',loc='top',fontsize=15)
            ax[gpos].yaxis.set_major_locator(plt.MultipleLocator(2))
            ax[gpos].yaxis.set_minor_locator(plt.NullLocator())
            ax[gpos].xaxis.set_major_locator(plt.MultipleLocator(4))
            ax[gpos].xaxis.set_minor_locator(plt.NullLocator())
        #ax[gpos].set_xlim(-5,0)        
        #ax[gpos].set_ylim(-2, 2)  
        #ax[gpos].set_zlim(-20, 20)         

        #cs=ax[gpos].contour(np.array(X),np.array(Y),np.transpose(np.array(Z)), levels=10,colors="black",linewidths=0.5,linestyles=["solid","dashed"])
    #ax[len(grids)-1].set_xlabel('$x$',loc='right')
    #ax[0].set_ylabel('$k^2\\;[\\mathrm{GeV^2}]$',rotation='vertical',loc='top')
    #fig.set_figheight(2.75)
    #fig.set_figwidth(3*len(grids))
    #fig.colorbar(surf)
    if savefile!="":
        fig.savefig(savefile)
    else:
        plt.show()

if __name__=="__main__":
    main()
