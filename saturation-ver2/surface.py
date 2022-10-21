#! /usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import numpy as np

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
#fig, ax = plt.subplots()

#ax.set(yscale="log",xscale="log")


#with open("./KSgluon-2.0/grids/KSnonlinear.dat","r") as fi:
#with open("../Run/GBW/Mass0.0-Qup650-Model0-Sud0/gluon-grid.dat","r") as fi:
with open("../Run/GBWS-Fix-S/Mass0.0-Qup650-Model22-Sud1/gluon-grid.dat","r") as fi:
#with open("../Run/BGK/Mass0.0-Qup650-Model1-Sud0/gluon-grid.dat","r") as fi:
#with open("../Run/BGKS-Fix-S/Mass0.0-Qup650-Model3-Sud1/gluon-grid.dat","r") as fi:
    X=[]
    Y=[]
    z=[]
    Z=[]
    x="n"
    counter=0
    posQ=49
    lineno=0
    prevQ2=0
    for i in fi:
        data=i.strip().split()
        #print(data)
        if len(data)==4:

            if float(prevQ2)>float(data[2]):
                lineno=0;
            lineno+=1;
            prevQ2=data[2]
            if(posQ!=(lineno-1)):
                continue
            else:
                print(data[2]);

        if((x=="n") or (data[0]!=x)):
            X.append(np.log10(np.exp(float(data[0]))))
            #X.append(np.exp(float(data[0])))
            if x!="n":
                #print(len(z))
                Z.append(z)
                z=[]
                counter+=1
            #z.append( float(data[len(data)-1]) )
            z.append( np.log10(abs(float(data[len(data)-1]))) )
            x=data[0];
        else:            
            z.append(np.log10( abs(float(data[len(data)-1])) ))
            #z.append( float(data[len(data)-1]) )
        if(counter==0):
            Y.append(np.log10(np.exp(float(data[1]))))
            #Y.append(np.exp(float(data[1])))
    Z.append(z)

print(np.array(X))
print(np.array(Y))
print(np.array(Z))
print(np.array(X).shape)
print(np.array(Y).shape)
print(np.array(Z).shape)
X, Y = np.meshgrid(X,Y )
#for i in range(len(Z)):
#    print(Z[i])
            


# Make data.
#X = np.arange(-5, 5, 0.25)
#Y = np.arange(-5, 5, 0.25)
#X, Y = np.meshgrid(X, Y)
#R = np.sqrt(X**2 + Y**2)
#Z = np.sin(R)

# Plot the surface.
#surf = ax.plot_surface(np.array(X),np.array(Y),np.transpose(np.array(Z)), cmap=cm.coolwarm,linewidth=0, antialiased=False)
#surf = ax.contourf(np.array(X),np.array(Y),np.transpose(np.array(Z)),cmap=cm.coolwarm)
surf = ax.plot_wireframe(np.array(X),np.array(Y),np.transpose(np.array(Z)), rcount=5, ccount=5)


# Customize the z axis.
#ax.set_zlim(-1.01, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
#ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
