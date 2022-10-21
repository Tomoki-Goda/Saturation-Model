#! /usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np
import getopt,sys



#fig, ax = plt.subplots()

def main():
    DIR="N/A"
    posQ2=0
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "s:g:q:", ["save=","grid=","Q2="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
    for opt, arg in opts:
        if opt in ["-s","--save"]:
            savefile=arg
        if opt in ["-g","--grid"]:
            DIR=arg
        if opt in ["-q","--Q2"]:
            print(arg)
            posQ2=int(arg)
    if DIR=="N/A":
        DIR=args[0]

    #DIR="../GKSgluon-1.1/grids/"
    #with open("./KSgluon-2.0/grids/KSnonlinear.dat","r") as fi:
    #with open("../Run/GBW/Mass0.0-Qup650-Model0-Sud0/gluon-grid.dat","r") as fi:
    #with open("../Run/GBWS-Fix-S/Mass0.0-Qup650-Model22-Sud1/gluon-grid.dat","r") as fi:
    #with open(DIR+"KSnonlinear.dat","r") as fi:
    with open(DIR,"r") as fi:
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
                #X.append(np.log10(np.exp(float(data[0]))))
                X.append(np.exp(float(data[0])))
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
                #Y.append(np.log10(np.exp(float(data[1]))))
                Y.append(np.exp(float(data[1])))
        Z.append(z)
    if len(Q2)>0:
        Q2=[np.log10(np.exp(float(i))) for i in list(Q2)]
        Q2.sort()
        print(Q2[0],"  ", Q2[len(Q2)-1]);
    print(X[0],"  ", X[len(X)-1]);
    print(Y[0],"  ", Y[len(Y)-1]);
    X, Y = np.meshgrid(X,Y )

    #fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    # Plot the surface.
    #surf = ax.plot_surface(np.array(X),np.array(Y),np.transpose(np.array(Z)), cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #surf = ax.plot_wireframe(np.array(X),np.array(Y),np.transpose(np.array(Z)), rcount=5, ccount=5)

    fig, ax = plt.subplots()
    ax.set(yscale="log",xscale="log")

    surf = ax.contourf(np.array(X),np.array(Y),np.transpose(np.array(Z)),levels=10,cmap=cm.coolwarm)
    cs=ax.contour(np.array(X),np.array(Y),np.transpose(np.array(Z)), levels=10,colors="black",linewidths=0.5,linestyles=["solid","dashed"])
    # Customize the z axis.
    #ax.set_zlim(-1.01, 1.01)
    #ax.zaxis.set_major_locator(LinearLocator(10))
    # A StrMethodFormatter is used automatically
    #ax.zaxis.set_major_formatter('{x:.02f}')

    # Add a color bar which maps values to colors.
    bar=fig.colorbar(surf, shrink=0.5, aspect=10)
    bar.add_lines(cs)

    plt.show()

if __name__=="__main__":
    main()
