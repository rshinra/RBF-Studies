# -*- coding: utf-8 -*-
"""
Initial attempt at Radial Basis Solving
(Solving Differential Equations using Radial Basis functions)
Function to model: u(x) = e^-x
Differential equation: du/dx = -u | u(0) = 1
NOTE: Didn't work for low numbers of centers, weirdo betas (large, small)
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
beta = np.linspace(1,20,39) # beta is actually Beta squared
err = zeros(np.size(beta))
    # "Test Data points" (function evaluation)
xd = np.linspace(0,1,10)
  # Missing Parameter for RBF Level Ia
#coeff = np.ones(np.size(xc)) #Initialized previously, but created from lstsq
  # Boundary condition to "solve problem"
xbound = 0
# Functions--------------------------------------------------------------------
 # Function for evaluation: f(x) -> x (vector)
def f(x):
   return np.exp(x*x/2)
   
 # Function Derivative for "test Data" f'(x) -> 1
def f_prime(x):
    return x*np.exp(x*x/2)
    
 # Basis Evaluation
def basis(i,x): # RADIAL Basis i evaluated at x
    return ( exp(-beta2 * (x-xc[i])*(x-xc[i]) ) )
    
 # Differential of basis; Recall, A Differential of a SUM is a SUM of Differentials
def basis_prime(i,x):
    return( 2 * beta2 * (xc[i]-x) * exp( -beta2 * (x-xc[i])*(x-xc[i]) ) )
    
 # Plot Final Function Approximation
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
    
# Plot Particular Basis Function
def plot_rbf(p):
    y = np.zeros(np.size(xd))
    for i in range(np.size(xd)):
        y[i] = coeff[p]*basis(p,xd[i])
    plt.plot(xd,y,'r')
def calc_error():
    y = np.zeros(np.size(xd))
    for i in range(np.size(xd)):
        for j in range(np.size(xc)):
            y[i] = y[i]+coeff[j]*basis(j,xd[i])
    error = np.zeros(np.size(xd))
    for i in range(np.size(xd)):
        error[i] = (ans[i]-y[i])*(ans[i]-y[i])
    return(np.sum(error))    
#------------------------------------------------------------------------------
# Main-------------------------------------------------------------------------
ans = f(xd)
for l in range(np.size(beta)):    
    plt.clf()
    beta2 = beta[l]*beta[l]
    # Build Matrix A
    A = np.zeros((np.size(xd)+1,np.size(xc)))
    for j in range(np.size(xc)): # J is columns of A
        for i in range(np.size(xd)): # I is rows of A
        # A for PDE is Differential equation with Basis and Basis_prime (j) evaluated at xd i
        # Differential equation to solve is: u' = 1 | รป' = 1 | A = phi'(x) g = 1
            A[i,j]=basis_prime(j,xd[i])
    # Build Matrix (Vector) (adding one for initial/boundary condition)
    
    g=np.zeros(np.size(xd)+1)
    for i in range(np.size(xd)):
        g[i] = f_prime(xd[i])
    
    # Build out Matrix and vector for boundary equation
    for j in range(np.size(xc)):
        A[np.size(xd),j]=basis(j,xbound)
    g[np.size(xd)]=f(xbound)
    # Solve
    a = np.linalg.lstsq(A,g)
    coeff=a[0] # Copy to coefficients
    err[l] = calc_error()
    #for i in range(np.size(xc)):
    # plot_rbf(i)
plot(beta,err)    
#plot_f()
#plt.plot(xd,ans,'b')
