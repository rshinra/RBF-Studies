# -*- coding: utf-8 -*-
"""
Initial attempt at Radial Basis Solving
(Solving Differential Equations using Radial Basis functions)
Function to model:      e^px-1/e^p-1     p = a/c
Differential equation:  au' = cu'' | u(0) = 0 | u(1) = 1
David Norman
"""
# Imports----------------------------------------------------------------------
import numpy as np
from math import exp
import matplotlib.pyplot as plt
#Global variables--------------------------------------------------------------
  # Input parameters for RBF Level Ia
    # Basis center Locations / Number of Bases
xc=np.linspace(0,1,50)
    # Width parameter (How Fat your basis functions are)
#beta = np.linspace(1,20,39)  
beta = 10
    # "Test Data points" (function evaluation)
xd = np.linspace(0,1,101)
  # Missing Parameter for RBF Level Ia
#coeff = np.ones(np.size(xc)) #Initialized previously, but created from lstsq
  # Boundary condition to "solve problem"
xbound0 = 0.
xbound1 = 1.
d1=20.
d2=6.
p=d1/d2
# Functions--------------------------------------------------------------------
 # Function for evaluation:  f(x) -> 
def f(x):
   return (np.exp(p*x) - 1)/(np.exp(p) - 1)
   
 # Function Derivative for "test Data" f'(x) -> 
def f_prime(x):
    return p*(np.exp(p*x))/(np.exp(p) - 1)

 # Function Second Derivative for "test Data" f''(x) -> 1
def f_double_prime(x):
    return p*p*(np.exp(p*x))/(np.exp(p)-1)

    
 # Basis Evaluation
def basis(b,x,cntr):                     # RADIAL Basis i evaluated at x
    b_2 = b*b
    xntr = x-cntr
    xntr_2 = xntr*xntr
    return exp(-b_2*xntr_2 ) 
    
 # Differential of basis; Recall, A Differential of a SUM is a SUM of Differentials
def basis_prime(b,x,cntr):
    b_2 = b*b
    xntr = x-cntr
    xntr_2 = xntr*xntr
    return -2*b_2*xntr*exp(-b_2*xntr_2)

 # Second differential of basis; Recall, A Differential of a SUM is a SUM of Differentials
def basis_double_prime(b,x,cntr):
    b_2 = b*b
    xntr = x-cntr
    xntr_2 = xntr*xntr
    return 2*b_2*exp(-b_2*xntr_2)*(2*b_2*xntr_2 - 1)
# TESTING A BUILD FUNCTION
def BuildAij(b,x,cntr):
    return(p*basis_prime(b,x,cntr)-basis_double_prime(b,x,cntr))

def BuildAMatrix(beta,cntrs):
    Amat = np.zeros((np.size(cntrs)+2,np.size(cntrs)))       # Why is there a 2 here?  (number of bc's)
    for j in range(np.size(cntrs)):     # J is columns of A
        for i in range(np.size(cntrs)):    # I is rows of A
            Amat[i,j]=BuildAij(beta,cntrs[i],cntrs[j])
    return Amat

 # Plot Final Function Approximation
def plot_f():
    y = np.zeros(np.size(xd))
    for i in range(np.size(xd)):
        for j in range(np.size(xc)):
            y[i] = y[i]+coeff[j]*basis(beta,xd[i],xc[j])
    plt.plot(xd,y,'k')
    yp = np.zeros(np.size(xc))
    for i in range(np.size(xc)):
        for j in range(np.size(xc)):
            yp[i] = yp[i]+coeff[j]*basis(beta,xc[i],xc[j])
    plt.plot(xc,yp,'ko')
    
# Plot Particular Basis Function
def plot_rbf(p):
    y = np.zeros(np.size(xd))
    for i in range(np.size(xd)):
        y[i] = coeff[p]*basis(beta,xd[i],xc[p])
    plt.plot(xd,y,'r')
#------------------------------------------------------------------------------
# Main-------------------------------------------------------------------------
plt.clf()
# Build Matrix A
A = BuildAMatrix(beta,xc)
# Build Matrix (Vector)  (adding one for initial/boundary condition)
g=np.zeros(np.size(xc)+2)                       # Same here with the "2"
# Build out Matrix and vector for boundary equations
for j in range(np.size(xc)):                #  Boundary Condition 0
    A[np.size(xc)-1,j]=basis(beta,xbound0,xc[j])
g[np.size(xc)-1]=f(xbound0)
for j in range(np.size(xc)):                #  Boundary Condition 1 (doesn't have to be at 1 and 0)
    A[np.size(xc),j]=basis(beta,xbound1,xc[j])
g[np.size(xc)]=f(xbound1)
# Solve
a = np.linalg.lstsq(A,g)
coeff=a[0]      # Copy to coefficients

ans = f(xd)
#for i in range(np.size(xc)):
#    plot_rbf(i)    
plot_f()
plt.plot(xd,ans,'b')
