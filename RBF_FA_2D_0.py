# -*- coding: utf-8 -*-
"""
Initial attempt at Radial Basis Interpolation
(Function Approximation using Radial Basis functions)
Function to model:        z = xy

David Norman
"""
# Imports
import numpy as np
from math import exp
import matplotlib.pyplot as plt

#Global variables
  # Input parameters for RBF Level Ia
    # Basis center Locations / Number of Bases
xc=np.linspace(0,1,10)
yc=np.linspace(0,1,10)
    # Width parameter (How Fat your basis functions are)
beta = 100                  # beta is actually Beta squared 
    # "Test Data" (function evaluation points)
xd = np.linspace(0,1,100)
yd = np.linspace(0,1,100)
  #Missing Parameter for RBF Level Ia is height parameter
coeff = np.zeros(np.size(xc)*np.size(yc)) # Placeholder

# Functions
 # Function for evaluation
#def f(x):
#   y=np.zeros(np.size(x))
#   for i in range(np.size(x)):
#       y[i]=xd[i]
#   return y
 # Basis Evaluation------------------------------------------
def basis(i,x,y):                     # Basis i evaluated at [x,y]
    return (exp(-beta*((x-xc[i])*(x-xc[i]))+(y-yc[i])*(y-yc[i]))
 # Plot Final Function Approximation-------------------------
def plot_f():
    z = np.zeros(np.size(xd)*np.size(yd))
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
plt.clf()
# Build Matrix A
A = np.zeros((np.size(xc),np.size(xc)))
for j in range(np.size(coeff)):     # J is columns of A
    for i in range(np.size(xc)):    # I is rows of A
    # A is Basis j evaluated at xd i
        A[i,j]=basis(j,xc[i])
# Build Matrix (Vector)
  #ALREADY BUILT, Y=X
#ANS = f(x)
# Solve
a = np.linalg.lstsq(A,np.transpose(xc))
coeff=a[0]

for i in range(np.size(xc)):
    plot_rbf(i)    
plot_f()
plt.plot(xd,xd,'b')
