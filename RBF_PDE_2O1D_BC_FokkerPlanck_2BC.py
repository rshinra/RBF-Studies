# -*- coding: utf-8 -*-
"""
Fokker Planck Attempt  POOR ATTEMPT, NO REAL BC's in FP equation
(Solving Differential Equations using Radial Basis functions)
Function to model:      r = -40x^3 + 10x | u = .455 e^[-10x^4+5x^2]
Differential equation:  (ru)' = u'' | u(-.5) = .85 | u(.5) = .85
David Norman
"""
# Imports----------------------------------------------------------------------
import numpy as np
from math import exp
import matplotlib.pyplot as plt
#Global variables--------------------------------------------------------------
  # Input parameters for RBF Level Ia
    # Basis center Locations / Number of Bases
xbound0 = -.5
xbound1 = .5
xc=np.linspace(xbound0,xbound1,50)
    # Width parameter (How Fat your basis functions are)
#beta = np.linspace(1,20,39)  
beta2 = 100
    # "Test Data points" (function evaluation)
xd = np.linspace(xbound0,xbound1,201)
  # Missing Parameter for RBF Level Ia
#coeff = np.ones(np.size(xc)) #Initialized previously, but created from lstsq
  # Boundary condition to "solve problem"
# Functions--------------------------------------------------------------------
 # Function for evaluation:  f(x) -> 
def f(x):
   return .455 * np.exp(-10*x*x*x*x + 5*x*x)
   
 # Function Derivative for "test Data" f'(x) -> 
def f_prime(x):
    return .455 * np.exp(-10*x*x*x*x + 5*x*x) * (-40*x*x*x + 10*x)

 # Function Second Derivative for "test Data" f''(x) -> 1
def f_double_prime(x):
    return .455 * np.exp(-10*x*x*x*x + 5*x*x) * (-120*x*x + 10 + (-40*x*x*x + 10*x)*(-40*x*x*x + 10*x))
def gamma(x):
    return -40*x*x*x + 10*x
def gamma_prime(x):
    return -120*x*x + 10
    
 # Basis Evaluation
def basis(i,x):                     # RADIAL Basis i evaluated at x
    return ( exp(-beta2 * (x-xc[i])*(x-xc[i]) ) )
    
 # Differential of basis; Recall, A Differential of a SUM is a SUM of Differentials
def basis_prime(i,x):
    return( 2 * beta2 * (xc[i]-x) * exp( -beta2 * (x-xc[i])*(x-xc[i]) ) )

 # Second differential of basis; Recall, A Differential of a SUM is a SUM of Differentials
def basis_double_prime(i,x):
    return( 2*beta2*exp( -beta2 * (x-xc[i])*(x-xc[i]) ) * (2*beta2*(x-xc[i])*(x-xc[i]) - 1) )

    
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
#------------------------------------------------------------------------------
# Main-------------------------------------------------------------------------
plt.clf()
# Build Matrix A
A = np.zeros((np.size(xc)+2,np.size(xc)))       # Why is there a 2 here?  (number of bc's)
for j in range(np.size(xc)):     # J is columns of A
    for i in range(np.size(xc)):    # I is rows of A
    # A for PDE is Differential equation with Basis, Basis_prime, and Basis Double Prime (j) evaluated at xd i
    # Differential equation to solve is: REPLACEu' = 1 | รป' = 1 | A = phi'(x) g = 1
        A[i,j]=basis_double_prime(j,xc[i])-basis_prime(j,xc[i])*gamma(xc[i])-basis(j,xc[i])*gamma_prime(xc[i])
# Build Matrix (Vector)  (adding one for initial/boundary condition)

g=np.zeros(np.size(xc)+2)                       # Same here with the "2"
#for i in range(np.size(xd)):
#    g[i] = f_double_prime(xd[i])

# Build out Matrix and vector for boundary equations
for j in range(np.size(xc)):                #  Boundary Condition 0
    A[np.size(xc)-1,j]=basis(j,xbound0)
g[np.size(xc)-1]=f(xbound0)
for j in range(np.size(xc)):                #  Boundary Condition 1 (doesn't have to be at 1 and 0)
    A[np.size(xc),j]=basis(j,xbound1)
g[np.size(xc)]=f(xbound1)

# Solve
a = np.linalg.lstsq(A,g)
coeff=a[0]      # Copy to coefficients

ans = f(xd)
#for i in range(np.size(xc)):
#    plot_rbf(i)    
plot_f()
plt.plot(xd,ans,'b')
