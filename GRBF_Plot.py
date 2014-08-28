# -*- coding: utf-8 -*-
"""
Created on Sat Jun 21 16:49:47 2014

@author: Rufus
"""

import numpy as np
from math import exp
import matplotlib.pyplot as plt
plt.clf()               # Clear the figure
centers = 5             # Number of Bases (centered at xi's)
xi=[.25,.33,.57,.75,.9] # Vector of centers, xi's dimension should be "centers"
lam=[3,5,2,7,7]         # lambda, multiplyer of basis, or ck in Meadian notation
beta = 100              # Beta, fatness factor
x = np.linspace(0,1,1000)   # x vector -> evaluation points

def f(x):               # f(x) add up all bases*height and then plot it
    y = np.linspace(0,1,np.size(x))
    for i in range(np.size(x)):
        y[i] = 0    
        for j in range(centers):
            y[i] = y[i]+lam[j]*(exp(-beta*abs(x[i]-xi[j])*abs(x[i]-xi[j])))
    plt.plot(x,y)
def rbf(x,p):           # Calc/plot a single basis function
    y = np.linspace(0,1,np.size(x))
    for i in range(np.size(x)):
        y[i] = 0    
        y[i] = y[i]+lam[p]*(exp(-beta*abs(x[i]-xi[p])*abs(x[i]-xi[p])))
    plt.plot(x,y)
# for each basis (i = 1 to centers) plot each basis, then plot the function
rbf(x,0)
rbf(x,1)
rbf(x,2)
rbf(x,3)
rbf(x,4)
f(x)
