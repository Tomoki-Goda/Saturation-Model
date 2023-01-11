#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


df=pd.read_table("./gbwscatter.txt",delimiter="\t")
df.columns=['kappa2','kt2','beta','Q2','x','mf2','val']
print(df)

fig,ax=plt.subplots(1,3,layout="constrained")

ax[0].scatter(df['kappa2'],df['kt2'],s=0.1)
ax[0].set(yscale="log", xscale="log")

ax[1].scatter(df['kappa2'],df['beta'],s=0.1)
ax[1].set( xscale="log")
#plt.show()
ax[2].scatter(df['kt2'],df['beta'],s=0.1)
ax[2].set( xscale="log")


plt.show()

