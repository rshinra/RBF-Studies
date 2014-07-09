# -*- coding: utf-8 -*-
"""
Created on Sat Jun 21 16:49:47 2014

@author: Rufus
"""

import numpy as np
from math import exp
import matplotlib.pyplot as plt
clf()
centers = 5
xi=[.25,.33,.57,.75,.9]
lam=[3,5,2,7,7]
beta = 100
x = np.linspace(0,1,1000)

def f(x):
    y = np.linspace(0,1,np.size(x))
    for i in range(np.size(x)):
        y[i] = 0    
        for j in range(centers):
            y[i] = y[i]+lam[j]*(exp(-beta*abs(x[i]-xi[j])*abs(x[i]-xi[j])))
    plt.plot(x,y)
def rbf(x,p):
    y = np.linspace(0,1,np.size(x))
    for i in range(np.size(x)):
        y[i] = 0    
        y[i] = y[i]+lam[p]*(exp(-beta*abs(x[i]-xi[p])*abs(x[i]-xi[p])))
    plt.plot(x,y)

rbf(x,0)
rbf(x,1)
rbf(x,2)
rbf(x,3)
rbf(x,4)
f(x)
