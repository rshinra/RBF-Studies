# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 17:42:06 2014

@author: Rufus
"""
import numpy as np                      #Numpy for arrays and such
from math import exp                    #uhh... e?
import matplotlib.pyplot as plt         #2D Plotting
from mpl_toolkits.mplot3d import axes3d #3D Plotting
from matplotlib import cm               #Color Maps

x = [.1,.3,.5,.7,.9]
y = [.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5]
X,Y = np.meshgrid(x,y)
# without transpose
for i in range (np.size(x)):
    for j in range(np.size(y)):
        print(i,j,X[j,i],Y[j,i])

# with transpose
X = np.transpose(X)
Y = np.transpose(Y)
for i in range (np.size(x)):
    for j in range(np.size(y)):
        print(i,j,X[i,j],Y[i,j])
x=np.linspace(-10,10,100)
y=np.linspace(-10,10,100)
X,Y = np.meshgrid(x,y)
X=np.transpose(X)
Y=np.transpose(Y)
Z=sqrt(X*X+Y*Y)
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X,Y,Z)