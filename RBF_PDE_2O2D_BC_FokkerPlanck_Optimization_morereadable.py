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
def Basis_2D(b,x,y,cx,cy):              # RADIAL Basis centered at (cx,cy) evaluated at (x,y)
    b_2 = b*b
    x_c = x - cx
    x_c_2 = x_c*x_c
    y_c = y - cy
    y_c_2 = y_c*y_c
    return ( exp(-b_2 * (x_c_2 + y_c_2) ))
def Basis_2D_x(b,x,y,cx,cy):              # RADIAL Basis centered at (cx,cy) evaluated at (x,y)
    b_2 = b*b
    x_c = x - cx
    x_c_2 = x_c*x_c
    y_c = y - cy
    y_c_2 = y_c*y_c
    return ( -2*b_2*x_c*exp(-b_2 * (x_c_2 + y_c_2) ))
def Basis_2D_y(b,x,y,cx,cy):              # RADIAL Basis centered at (cx,cy) evaluated at (x,y)
    b_2 = b*b
    x_c = x - cx
    x_c_2 = x_c*x_c
    y_c = y - cy
    y_c_2 = y_c*y_c
    return ( -2*b_2*y_c*exp(-b_2 * (x_c_2 + y_c_2) ))
def Basis_2D_xx(b,x,y,cx,cy):              # RADIAL Basis centered at (cx,cy) evaluated at (x,y)
    b_2 = b*b
    x_c = x - cx
    x_c_2 = x_c*x_c
    y_c = y - cy
    y_c_2 = y_c*y_c
    return ( 2*b_2*(2*b_2*x_c_2 - 1)*exp(-b_2 * (x_c_2 + y_c_2) ))
def Basis_2D_xy(b,x,y,cx,cy):              # RADIAL Basis centered at (cx,cy) evaluated at (x,y)
    b_2 = b*b
    x_c = x - cx
    x_c_2 = x_c*x_c
    y_c = y - cy
    y_c_2 = y_c*y_c
    return ( 4*b_2*x_c*y_c*exp(-b_2 * (x_c_2 + y_c_2) ))
def Basis_2D_yy(b,x,y,cx,cy):              # RADIAL Basis centered at (cx,cy) evaluated at (x,y)
    b_2 = b*b
    x_c = x - cx
    x_c_2 = x_c*x_c
    y_c = y - cy
    y_c_2 = y_c*y_c
    return ( 2*b_2*(2*b_2*y_c_2 - 1)*exp(-b_2 * (x_c_2 + y_c_2) ))
def f(x,y):
    return(.00831521 * np.exp(-0.5*(-x*x + 0.05*x*x*x*x + y*y)))
def g1(x,y):
    return y
def g1_x(x,y):
    return 0
def g1_y(x,y):
    return 1
def g2(x,y):
    xi = .2
    w = 1
    gam = .1
    return - 2*xi*w*y + w*w*x - w*w*gam*x*x*x
def g2_x(x,y):
    xi = .2
    w = 1
    gam = .1
    return w*w - 3*w*w*gam*x*x
def g2_y(x,y):
    xi = .2
    w = 1
    gam = .1
    return - 2*xi*w
# Matrixes for the problem
def Build_A(b,x,y,xc,yc):               # Matrix H (Ax=b, Hc=g)
    sx = np.size(x,0)
    sy = np.size(y,1)
    sxy = sx*sy
    sxc = np.size(xc,0)
    syc = np.size(yc,1)
    sxyc = sxc*syc
    A = np.zeros(sxy+1,sxy+1)       # Why is there a 1 here?  (Integral Boundary Condition)
    for i in range(sx):     
        for j in range(sy):    # x[i](maxy) + y[j] is rows of A
            for k in range(sxc):
                for l in range(syc): # xc[i](maxy) + yc[j] is column of A
                    Ai = i*(sy)+j
                    Aj = k*(syc)+l      
                    A[Ai,Aj]=Aij(b,x[i,j],y[i,j],xc[k,l],yc[k,l])
    for i in range(sxy):
        A[i,sxyc]=1.         # for Lagrange Multiplier
    return A
def Build_b(l):
    b = np.zeros(l+1)
    return b
def Aij(b,x,y,xc,yc):
    D = .4
    return D*Basis_2D_yy(b,x,y,xc,yc)-Basis_2D(b,x,y,xc,yc)*g1_x(x,y)-g1(x,y)*Basis_2D_x(b,x,y,xc,yc)-g2_y(x,y)*Basis_2D(b,x,y,xc,yc)-g2(x,y)*Basis_2D_y(b,x,y,xc,yc)
def u_hat_x(b,x,y,xc,yc,c):
    sxc = np.size(xc,0)
    syc = np.size(yc,1)
    z = 0
    for i in range(sxc):
        for j in range(syc):
            z = z + c[i,j]*Basis_2D(b, x, y, xc[i,j],yc[i,j])
    return z
def u_hat(b,x,y,xc,yc,c):
    sx = np.size(x,0)
    sy = np.size(y,1)
    z = np.zeros((sx,sy))
    for i in range(sx):               #Columns of X, and Rows of Y
        for j in range(sy):           #For Each Z point (i,j):
            z[i,j] = u_hat_x(b,x[i,j],y[i,j],xc,yc,c)
    return z
def error(b,x,y,xc,yc,c):
    err = 0
    z = f(x,y)
    zapx = u_hat(b,x,y,xc,yc,c)
    errv = np.abs(z-zapx)
    err = np.sum(errv)    
    return err
def errfinal(b,x,y,xc,yc,c):
    err = 0
    z = f(x,y)
    f1 = plt.figure()
    ax = f1.gca(projection='3d')
    ax.plot_surface(x,y,z, rstride=8, cstride=8, alpha=0.3)
    zapx = u_hat(b,x,y,xc,yc,c)
    f2 = plt.figure()
    ax = f2.gca(projection='3d')
    ax.plot_surface(x,y,zapx, rstride=8, cstride=8, alpha=0.3)
    errv = np.abs(z-zapx)
    f3 = plt.figure()
    ax = f3.gca(projection='3d')
    ax.plot_surface(x,y,errv, rstride=8, cstride=8, alpha=0.3)
    err = np.sum(errv)    
    return err
def plot_f(b,x,y,xc,yc,c_m):               # f(x) add up all bases*height and then plot it
    z = f(x,y)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(x,y,z, rstride=8, cstride=8, alpha=0.3)
    cset = ax.contour(x,y,z, zdir='z', offset = -.5, cmap=cm.coolwarm)
    cset = ax.contour(x,y,z, zdir='x', offset = -6, cmap=cm.coolwarm)
    cset = ax.contour(x,y,z, zdir='y', offset = 6, cmap=cm.coolwarm)
    ax.set_xlim(-6,6)
    ax.set_ylim(-6,6)
    ax.set_zlim(-.5,.5)
    plt.show()
def plot_f_approx(ax,b,x,y,xc,yc,c):               # f(x) add up all bases*height and then plot it
    z = u_hat(b,x,y,xc,yc,c)
    ax.plot_surface(x,y,z, rstride=8, cstride=8, alpha=0.3)
    cset = ax.contour(x,y,z, zdir='z', offset = -.1, cmap=cm.coolwarm)
    cset = ax.contour(x,y,z, zdir='x', offset = -6, cmap=cm.coolwarm)
    cset = ax.contour(x,y,z, zdir='y', offset = 6, cmap=cm.coolwarm)
    ax.set_xlim(-6,6)
    ax.set_ylim(-6,6)
    ax.set_zlim(-.1,.1)
    plt.show()
def plot_f_error(ax,b,x,y,xc,yc,c):
    z = f(x,y)
    za = u_hat(b,x,y,xc,yc,c)
    ax.plot_surface(x,y,z-za, rstride=8, cstride=8, alpha=0.3)
    cset = ax.contour(x,y,z-za, zdir='z', offset = -.1, cmap=cm.coolwarm)
    cset = ax.contour(x,y,z-za, zdir='x', offset = -6, cmap=cm.coolwarm)
    cset = ax.contour(x,y,z-za, zdir='y', offset = 6, cmap=cm.coolwarm)
    ax.set_xlim(-6,6)
    ax.set_ylim(-6,6)
    ax.set_zlim(-.1,.1)
    plt.show()
# beta = .69
Zchk = 0
Zapxchk = 0
errvchk = 0
nbx = 10
nby = 10
nbasis = nbx*nby
nberx = 20
nbery = 20
nEEP = mberx*nbery
nbex = 201
nbey = 201
nEP = nbex*nbey
xbound0 = -2
xbound1 = 2
ybound0 = -2
ybound1 = 2
xbounde0 = -5
xbounde1 = 5
ybounde0 = -5
ybounde1 = 5
xc=np.linspace(xbound0,xbound1,nbx) # Vector of centers, xi's dimension should be "centers"
yc=np.linspace(ybound0,ybound1,nby) # Vector of centers, yi's dimension should be "centers"
Xc,Yc = np.meshgrid(xc,yc,indexing = 'ij')        # Please note: X is COLUMNS, Y is ROWS
xr = np.linspace(xbound0,xbound1,nberx)
yr = np.linspace(ybound0,ybound1,nbery)
Xr,Yr = np.meshgrid(xr,yr,indexing ='ij')
xe=np.linspace(xbounde0,xbounde1,nbex)
ye=np.linspace(ybounde0,ybounde1,nbey)
Xe,Ye = np.meshgrid(xe,ye,indexing ='ij')
betamin = .3
betamax = 32
betaspace = np.arange(betamin,betamax,.1)
plt.clf
minerr=1e9
lam = np.zeros(np.size(betaspace))
SumC = np.zeros(np.size(betaspace))
Errr = np.zeros(np.size(betaspace))
g = Build_b(nbx*nbx)
g[nbasis]=1
for i in range(np.size(betaspace)):
    beta = betaspace[i]    
    H = Build_A(beta,Xc,Yc,Xc,Yc)
#    g = Build_b(np.size(Xc[0])*np.size(Yc[1]))   (build every timestep if not constant length)
    # Boundary equations
    for j in range(nbasis):                #  Boundary Condition 0, integral, i.e. all ck must add up to that
      H[nbasis,j]=np.pi/beta/beta
    # Solve
    a = np.linalg.lstsq(H,g)
    ck=a[0]      # Copy to coefficients
    if (np.abs(ck[nbasis]) > .000001):
      lam[i] = ck[nbasis]
    else:
      lam[i] = 0
    SumC[i] = sum(H[nbasis]*ck)
    ckec = np.resize(ck[0:nbasis-1],(nbasis))
    Errr[i] = error(beta,Xr,Yr,Xc,Yc,ckec)
    if(Errr[i] < minerr):
        betafinal = beta
        ckfinal = ck
        lamfinal = lam[i]
        SCfinal = SumC[i]
        minerr = Errr[i]
        
#for j in range(np.size(Xc[0])*np.size(Yc[1])):                #  Boundary Condition 0, integral, i.e. all ck must add up to that
#    H[np.size(Xc[0])*np.size(Yc[1]),j]=np.pi/beta/beta
#g[np.size(Xc[0])*np.size(Yc[1])]=1
#a = np.linalg.lstsq(H,g)
#ck=a[0]      # Copy to coefficients
lam=ckfinal[nbasis]
print(betafinal,lam,minerr)
ck = np.resize(ckfinal[0:nbasis-1],(nbx,nby))
fig = plt.figure(figsize=plt.figaspect(0.5))
ax1 = fig.add_subplot(1,2,1,projection='3d')
plot_f_approx(ax1,betafinal,Xe,Ye,Xc,Yc,ck)
#ax1 = plt.subplot(221)
#ax2 = plt.subplot(223)
#ax3 = plt.subplot(122)
ax2 = fig.add_subplot(1,2,2, projection='3d')
plot_f_error(ax2,betafinal,Xe,Ye,Xc,Yc,ck)
#fig3 = plt.figure
#plot_f(fig3,betafinal,xd,xc,ckfinal)
#ax1.plot(xd,ans,'b')
fig2 = plt.figure()
plt.plot(betaspace,Errr,'b')
plt.ylabel("ErrorSum")
plt.xlabel("Fatness Factor (Width Parameter)")
plt.ylim(0,2.5)