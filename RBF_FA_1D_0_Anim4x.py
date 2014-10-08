# -*- coding: utf-8 -*-
"""
Initial attempt at Radial Basis Interpolation
(Function Approximation using Radial Basis functions)
Function to model:        y = x

David Norman
"""
# Imports
import time
import numpy as np
from math import exp
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#Global variables
  # Input parameters for RBF Level ~IIb
    # Width parameter (How Fat your basis functions are)
beta = 100                  # beta is actually Beta squared 
    # "Test Data" (function evaluation points)
xd = np.linspace(0,1,1000)
  #Missing Parameter for RBF Level ~IIb is centers, animation of centers will show addition of bases helps fit
#xc=np.linspace(0,1,26)
  #Missing Parameter for RBF Level ~IIb is height parameter
#coeff = np.zeros(np.size(xc)) # Placeholder
# For the animation, we'll go from 1 to size(xc), adding one basis, and plotting the result, attempting to make an animation of this:

# Functions
 # Function for evaluation
#def f(x):
#   y=np.zeros(np.size(x))
#   for i in range(np.size(x)):
#       y[i]=xd[i]
#   return y
 # Basis Evaluation------------------------------------------
def basis(i,x):                     # Basis i evaluated at x
    return (exp(-beta*(x-xc[i])*(x-xc[i])))
 # Plot Final Function Approximation-------------------------
def plot_f():
    y = np.zeros(np.size(xd))
    for i in range(np.size(xd)):
        for j in range(np.size(xc)):
            y[i] = y[i]+coeff[j]*basis(j,xd[i])
    plt.plot(xd,y,'k')
    yp = np.zeros(np.size(xc))
    for i in range(np.size(xc)):
        for j in range(np.size(xc)):
            yp[i] = yp[i]+coeff[j]*basis(j,xc[i])
    plt.plot(xc,yp,'ko')
# Plot Particular Basis Function-------------------
def plot_rbf(p):
    y = np.zeros(np.size(xd))
    for i in range(np.size(xd)):
        y[i] = coeff[p]*basis(p,xd[i])
    plt.plot(xd,y,'r')

# Main
for b in range(26):
    fig1=plt.figure()
    plt.clf()
    xc=np.linspace(0,1,b+1)
    coeff = np.zeros(np.size(xc))
    A = np.zeros((np.size(xc),np.size(xc)))
    for j in range(b+1):     # J is columns of A
        for i in range(np.size(xc)):    # I is rows of A
        # A is Basis j evaluated at xd i
            A[i,j]=basis(j,xc[i])
    # Build Matrix (Vector)
      #ALREADY BUILT, Y=X
    #ANS = f(x)
    # Solve
    a = np.linalg.lstsq(A,np.transpose(xc))
    coeff=a[0]
    
    #for i in range(b+1):
        #plot_rbf(i)    
    plot_f()
    plt.plot(xd,xd,'b')
    plt.draw
    plt.pause(.0001)
    time.sleep(0.5)

