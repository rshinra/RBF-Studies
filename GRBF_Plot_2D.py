# -*- coding: utf-8 -*-
"""
Created on Sat Jun 21 16:49:47 2014

@author: Rufus
"""

import numpy as np                      #Numpy for arrays and such
from math import exp                    #uhh... e?
import matplotlib.pyplot as plt         #2D Plotting
from mpl_toolkits.mplot3d import axes3d #3D Plotting
from matplotlib import cm               #Color Maps
beta = 5
xc=np.linspace(-1,1,11) # Vector of centers, xi's dimension should be "centers"
yc=np.linspace(-1,1,11) # Vector of centers, yi's dimension should be "centers"
Xc,Yc = np.meshgrid(xc,yc)        # Please note: X is COLUMNS, Y is ROWS
xe=np.linspace(-2,2,151)
ye=np.linspace(-2,2,151)
Xe,Ye = np.meshgrid(xe,ye)
Ck=np.ones(np.size(xc) * np.size(yc))         # lambda, multiplyer of basis, or ck in Meadian notation
Ck=np.resize(Ck,(np.size(xc),np.size(yc)))             # Beta, fatness factor
def Basis_2D(b,x,y,cx,cy):              # RADIAL Basis centered at (cx,cy) evaluated at (x,y)
    b_2 = b*b
    x_c = x - cx
    x_c_2 = x_c*x_c
    y_c = y - cy
    y_c_2 = y_c*y_c
    return ( exp(-b_2 * (x_c_2 + y_c_2) ))  
def U_x(b,x,y,Xc,Yc,c):
    z = 0
    for i in range(np.size(Xc,1)):
        for j in range(np.size(Yc,0)):
            z = z + c[j,i]*Basis_2D(b, x, y, Xc[j,i],Yc[j,i])
    return z
def u_hat(b,x,y,xc,yc,c):
    z = np.zeros((np.size(x,1),np.size(y,0)))
    for i in range(np.size(x,1)):               #Columns of X, and Rows of Y
        for j in range(np.size(y,0)):           #For Each Z point (i,j):
            z[j,i] = U_x(b,x[j,i],y[j,i],xc,yc,c)
    return z
#def u_hat(b,x,y,xc,yc,c):
#    z = np.zeros((np.size(x,1),np.size(y,0)))
#    for i in range(np.size(x,1)):               #Columns of X, and Rows of Y
#        for j in range(np.size(y,0)):           #For Each Z point (i,j):
#            for k in range(np.size(xc,1)):
#                for l in range(np.size(yc,0)):
#                    z[j,i] = z[j,i] + z[j,i] + c[l,k]*Basis_2D(b, x[j,i], y[j,i], xc[l,k],yc[l,k])
#    return z
def plot_f(b,x,y,xc,yc,c_m):               # f(x) add up all bases*height and then plot it
    z = u_hat(b,x,y,xc,yc,c_m)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(x,y,z, rstride=8, cstride=8, alpha=0.3)
def error(b,x,y,xc,yc,c):
    err = 0
    z = f(x,y)
    errv = np.abs(z-u_hat(b,x,y,xc,yc,c))
    err = np.sum(errv)    
    return err
 # Plot Final Function Approximation
#def plot_f(ax,b,x,c_v,c_k):
#    ax.plot(x,u_hat(b,x,c_v,c_k),'k')
#    ax.plot(c_v,u_hat(b,c_v,c_v,c_k),'ko')
def rbf(lam):           # Calc/plot a single basis function
    x = np.linspace(-1,1,100)
    y = np.linspace(-1,1,100)
    Y,X = np.meshgrid(x,y)
    Z = np.zeros((np.size(x),np.size(y)))    
    for i in range(np.size(x)):
        for j in range(np.size(y)):
            Z[i,j] = Basis_2D(beta,X[i,j],Y[i,j],0,0)
            print(round(X[i,j],2),round(Y[i,j],2),round(Z[i,j],2))
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(X,Y,Z)
plot_f(beta,Xc,Yc,Xc,Yc,Ck)
#plot_f(beta,Xe,Ye,Xc,Yc,Ck)
#rbf(1)
#x = np.linspace(0,1,100)
#y = np.linspace(0,1,100)
#X,Y=np.meshgrid(x,y)
#f(X,Y)
