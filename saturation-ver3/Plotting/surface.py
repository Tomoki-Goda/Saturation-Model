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
    savefile="./grid-plot.jpg"
    dipole=False
    contour=False
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "dmhs:g:q:c", ['dipole','multiple','help',"save=","grid=","Q2=","contour"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
    for opt, arg in opts:
        if opt in ["-s","--save"]:
            savefile=arg
        if opt in ["-g","--grid"]:
            grids=[arg]
        if opt in ["-q","--Q2"]:
            print(arg)
            posQ2=int(arg)
        if opt in ['-d','--dipole']:
            dipole=True
        if opt in ['-h', '--hellp']:
            print("-s <save file>\n-g <grid>\n-q <position of Q2 in the list of Q2 in the grid...> ")
            print("-g flag is no longer needed simply write grid file at the end")
            print("-m to say multiple grid files, and list grids files after options. ")
            sys.exit()
        if opt in ['-m','--multiple']:
            grids=args;
        if opt in ['-c','--contour']:
            contour=True
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
    	fig, ax =plt.subplots(1,len(grids), sharex=True,sharey=True,layout='constrained',subplot_kw={'projection':'3d'})
    
    
    fig.set_constrained_layout_pads(hspace=0,wspace=0,h_pad=0.025,w_pad=0.025)
    if len(grids)==1:
        ax=[ax]
    for gpos in range(len(grids)):
        grid=grids[gpos]
        with open(grid,"r") as fi:
        #with open("../Run/BGK/Mass0.0-Qup650-Model1-Sud0/gluon-grid.dat","r") as fi:
        #with open("../Run/BGKS-Fix-S/Mass0.0-Qup650-Model3-Sud1/gluon-grid.dat","r") as fi:
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
                #print(data)
                #print(data)
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

                if((x=="n") or (data[0]!=x)):
                    X.append(np.log10(float(data[0])))
                    #X.append(np.exp(float(data[0])))
                    #X.append(float(data[0]))
                    if x!="n":
                        #print(len(z))
                        Z.append(z)
                        z=[]
                        counter+=1
                    z.append( float(data[len(data)-1]) )
                    #z.append( np.log10(abs(float(data[len(data)-1]))) )
                    x=data[0];
                else:            
                    #z.append(np.log10( abs(float(data[len(data)-1])) ))
                    z.append( float(data[len(data)-1]) )
                if(counter==0):
                    Y.append(np.log10(float(data[1])))
                    #Y.append(np.exp(float(data[1])))
                    #Y.append(float(data[1]))

            Z.append(z)
        if len(Q2)>0:
            Q2=[np.log10(np.exp(float(i))) for i in list(Q2)]
            Q2.sort()
            print(Q2[0],"  ", Q2[len(Q2)-1]);
        print(X[0],"  ", X[len(X)-1]);
        print(Y[0],"  ", Y[len(Y)-1]);
        X, Y = np.meshgrid(X,Y )
        
        #ax[gpos].set(yscale="log",xscale="log",xlim=[1.0e-7,1.0e-2],ylim=[0.05,5.0e+2])
        if contour:
            surf = ax[gpos].contourf(np.array(X),np.array(Y),np.transpose(np.array(Z)),levels=10,cmap=cm.coolwarm)
        else:
            surf = ax[gpos].plot_wireframe(np.array(X),np.array(Y),np.transpose(np.array(Z)), rcount=6, ccount=0,color="r",ls="-." )
            surf = ax[gpos].plot_wireframe(np.array(X),np.array(Y),np.transpose(np.array(Z)), rcount=0, ccount=6,color="b",ls="-")
            if dipole:
                ax[gpos].view_init(28, -28)
                ax[gpos].set_ylabel('$\log_{10}(r^2\\;[\\mathrm{GeV^2}])$',rotation='vertical',loc='top')
            else:
                ax[gpos].view_init(30, 25)
                ax[gpos].set_ylabel('$\log_{10}(k^2\\;[\\mathrm{GeV^2}])$',rotation='vertical',loc='top')
            ax[gpos].set_xlabel('$\log_{10}x$',loc='right')
            ax[gpos].set_zlim3d(-0.01, 0.001)        
        #cs=ax[gpos].contour(np.array(X),np.array(Y),np.transpose(np.array(Z)), levels=10,colors="black",linewidths=0.5,linestyles=["solid","dashed"])
    #ax[len(grids)-1].set_xlabel('$x$',loc='right')
    #ax[0].set_ylabel('$k^2\\;[\\mathrm{GeV^2}]$',rotation='vertical',loc='top')
    #fig.set_figheight(2.75)
    #fig.set_figwidth(3*len(grids))
    #fig.set_figheight(10)
    fig.set_figwidth(10*len(grids))
    
    fig.savefig(savefile)
    plt.show()

if __name__=="__main__":
    main()
