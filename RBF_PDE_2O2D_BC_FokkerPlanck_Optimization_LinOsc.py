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
#def f(x,y):
#    return(.00831521 * np.exp(-0.5*(-x*x + 0.05*x*x*x*x + y*y)))
def g1(x,y):
    return y
def g1_x(x,y):
    return 0
def g1_y(x,y):
    return 1
def g2(x,y):
    xi = .05
    w = 1
    gam = 0
    return - 2*xi*w*y - w*w*x - w*w*gam*x*x*x
def g2_x(x,y):
    xi = .05
    w = 1
    gam = 0
    return -w*w - 3*w*w*gam*x*x
def g2_y(x,y):
    xi = .05
    w = 1
    gam = 0
    return - 2*xi*w
# Matrixes for the problem
def Build_A(b,x,y,xc,yc):               # Matrix H (Ax=b, Hc=g)
    A = np.zeros((np.size(x,0)*np.size(y,1)+1,np.size(xc,0)*np.size(yc,1)+1))       # Why is there a 1 here?  (Integral Boundary Condition)
    for i in range(np.size(x,0)):     
        for j in range(np.size(y,1)):    # x[i](maxy) + y[j] is rows of A
            for k in range(np.size(xc,0)):
                for l in range(np.size(yc,1)): # xc[i](maxy) + yc[j] is column of A
                    Ai = i*(np.size(y,1))+j
                    Aj = k*(np.size(yc,1))+l      
                    A[Ai,Aj]=Aij(b,x[i,j],y[i,j],xc[k,l],yc[k,l])
    for i in range(np.size(x,0)*np.size(y,1)):
        A[i,np.size(xc,0)*np.size(yc,1)]=1.         # for Lagrange Multiplier
    return A
def Build_b(l):
    b = np.zeros(l+1)
    return b
def Aij(b,x,y,xc,yc):
    D = .1
    return D*Basis_2D_yy(b,x,y,xc,yc)-Basis_2D(b,x,y,xc,yc)*g1_x(x,y)-g1(x,y)*Basis_2D_x(b,x,y,xc,yc)-g2_y(x,y)*Basis_2D(b,x,y,xc,yc)-g2(x,y)*Basis_2D_y(b,x,y,xc,yc)
def u_hat_x(b,x,y,Xc,Yc,c):
    z = 0
    for i in range(np.size(Xc,0)):
        for j in range(np.size(Yc,1)):
            z = z + c[i,j]*Basis_2D(b, x, y, Xc[i,j],Yc[i,j])
    return z
def u_hat(b,x,y,xc,yc,c):
    z = np.zeros((np.size(x,0),np.size(y,1)))
    for i in range(np.size(x,0)):               #Columns of X, and Rows of Y
        for j in range(np.size(y,1)):           #For Each Z point (i,j):
            z[i,j] = u_hat_x(b,x[i,j],y[i,j],xc,yc,c)
    return z
#def error(b,x,y,xc,yc,c):
#    err = 0
#    z = f(x,y)
#    errv = (z-u_hat(b,x,y,xc,yc,c))
#    errv = errv*errv
#    err = np.sqrt(np.sum(errv))    
#    return err
def merit(b,x,y,xc,yc,c,l):
    z = u_hat(b,x,y,xc,yc,c)
    basement = np.sum(np.abs((z - np.abs(z))/2))
    attic = np.sum(((z-1) + np.abs(z-1))/2)
    return np.abs(l) + basement + attic
#def plot_f(ax,x,y):               # f(x) add up all bases*height and then plot it
#    z = f(x,y)
#    ax.plot_surface(x,y,z, rstride=8, cstride=8, alpha=0.3)
#    cset = ax.contour(x,y,z, zdir='z', offset = -np.amax(np.abs(z)), cmap=cm.coolwarm)
#    cset = ax.contour(x,y,z, zdir='x', offset = -np.amin(x), cmap=cm.coolwarm)
#    cset = ax.contour(x,y,z, zdir='y', offset = np.amax(y), cmap=cm.coolwarm)
#    ax.set_xlim(np.amin(x),np.amax(x))
#    ax.set_ylim(np.amin(y),np.amax(y))
#    ax.set_zlim(-round(np.amax(np.abs(z)),1),round(np.amax(np.abs(z)),1))
#    plt.show()
def plot_f_approx(ax,b,x,y,xc,yc,c):               # f(x) add up all bases*height and then plot it
    z = u_hat(b,x,y,xc,yc,c)
    ax.plot_surface(x,y,z, rstride=8, cstride=8, alpha=0.3)
    cset = ax.contour(x,y,z, zdir='z', offset = -np.amax(np.abs(z)), cmap=cm.coolwarm)
    cset = ax.contour(x,y,z, zdir='x', offset = np.amin(x), cmap=cm.coolwarm)
    cset = ax.contour(x,y,z, zdir='y', offset = np.amax(y), cmap=cm.coolwarm)
    ax.set_xlim(np.amin(x),np.amax(x))
    ax.set_ylim(np.amin(y),np.amax(y))
    ax.set_zlim(-round(np.amax(np.abs(z)),1),round(np.amax(np.abs(z)),1))
    plt.show()
#def plot_f_error(ax,b,x,y,xc,yc,c):
#    z = f(x,y)
#    za = u_hat(b,x,y,xc,yc,c)
#    z = z-za
#    ax.plot_surface(x,y,z, rstride=8, cstride=8, alpha=0.3)
#    cset = ax.contour(x,y,z, zdir='z', offset = -np.amax(np.abs(z)), cmap=cm.coolwarm)
#    cset = ax.contour(x,y,z, zdir='x', offset = np.amin(x), cmap=cm.coolwarm)
#    cset = ax.contour(x,y,z, zdir='y', offset = np.amax(y), cmap=cm.coolwarm)
#    ax.set_xlim(np.amin(x),np.amax(x))
#    ax.set_ylim(np.amin(y),np.amax(y))
#    ax.set_zlim(-.1,.1)
#    plt.show()
nbx = 11
nby = 11
xbound0 = -5
xbound1 = 5
ybound0 = -5
ybound1 = 5
xc=np.linspace(xbound0,xbound1,nbx) # Vector of centers, xi's dimension should be "centers"
yc=np.linspace(ybound0,ybound1,nby) # Vector of centers, yi's dimension should be "centers"
Xc,Yc = np.meshgrid(xc,yc,indexing = 'ij')        # Please note: X is COLUMNS, Y is ROWS
xr = np.linspace(-10,10,21)
yr = np.linspace(-10,10,21)
Xr,Yr = np.meshgrid(xr,yr,indexing ='ij')
xe=np.linspace(-10,10,201)
ye=np.linspace(-10,10,201)
Xe,Ye = np.meshgrid(xe,ye,indexing ='ij')
g = Build_b(np.size(Xc,0)*np.size(Yc,1))
g[np.size(Xc,0)*np.size(Yc,1)]=1
#Coarse
betamin = 1.135
betamax = 1.138
betaspace = np.arange(betamin,betamax,.00005)
lam = np.zeros(np.size(betaspace))
SumC = np.zeros(np.size(betaspace))
merits = np.zeros(np.size(betaspace))
#Errr = np.zeros(np.size(betaspace))
minmerit=1e9
for i in range(np.size(betaspace)):
    beta = betaspace[i]    
    H = Build_A(beta,Xc,Yc,Xc,Yc)
#    g = Build_b(np.size(Xc[0])*np.size(Yc[1]))   (build every timestep if not constant length)
    # Boundary equations
    for j in range(np.size(Xc,0)*np.size(Yc,1)):                #  Boundary Condition 0, integral, i.e. all ck must add up to that
      H[np.size(Xc,0)*np.size(Yc,1,),j]=np.pi/beta/beta
    # Solve
#    a = np.linalg.lstsq(H,g)
#    ck=a[0]      # Copy to coefficients
    ck = np.linalg.solve(H,g)
    if (np.abs(ck[np.size(Xc,0)*np.size(Yc,1)]) > .000001):
      lam[i] = ck[np.size(Xc,0)*np.size(Yc,1)]
    else:
      lam[i] = 0
    SumC[i] = sum(H[np.size(Xc,0)*np.size(Yc,1)]*ck)
    ckec = np.resize(ck[0:np.size(Xc,0)*np.size(Yc,1)-1],(np.size(Xc,0),np.size(Yc,1)))
    merits[i] = merit(beta,Xr,Yr,Xc,Yc,ckec,ck[np.size(Xc,0)*np.size(Yc,1)])
#    Errr[i] = error(beta,Xr,Yr,Xc,Yc,ckec)
    if(merits[i] < minmerit):
        betafinal = beta
        ckfinal = ck
        lamfinal = lam[i]
        SCfinal = SumC[i]
        minmerit = merits[i]
Bx = np.zeros((np.size(Xc,0)*np.size(Yc,1),1))
By = np.zeros((np.size(Xc,0)*np.size(Yc,1),1))
for i in range(np.size(Xc,0)):
    for j in range(np.size(Yc,1)):
        Bi = i*(np.size(Yc,1))+j
        Bx[Bi,0] = Xc[i,j]
        By[Bi,0] = Yc[i,j]
Bck = ck[0:np.size(Xc,0)*np.size(Yc,1)-1]
ck = np.resize(ckfinal[0:np.size(Xc,0)*np.size(Yc,1)-1],(np.size(Xc,0),np.size(Yc,1)))
print(betafinal,minmerit,np.amax(u_hat(betafinal,Xr,Yr,Xc,Yc,ck)))
fig = plt.figure(figsize=plt.figaspect(0.5))
ax1 = fig.gca(projection='3d')
#ax1 = fig.add_subplot(1,2,1,projection='3d')
plot_f_approx(ax1,betafinal,Xe,Ye,Xc,Yc,ck)
#ax1 = plt.subplot(221)
#ax2 = plt.subplot(223)
#ax3 = plt.subplot(122)
#ax2 = fig.add_subplot(1,2,2, projection='3d')
#plot_f_error(ax2,betafinal,Xe,Ye,Xc,Yc,ck)
#plot_f(ax2,Xe,Ye)
#fig3 = plt.figure
#plot_f(fig3,betafinal,xd,xc,ckfinal)
#ax1.plot(xd,ans,'b')
#Plot error over betaspace
fig2 = plt.figure()
plt.plot(betaspace,merits,'b')
plt.ylabel("Merit")
plt.xlabel("Fatness Factor (Width Parameter)")
plt.ylim(0,np.amax(merits))