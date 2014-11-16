# -*- coding: utf-8 -*-
"""
Solving Convection diffusion equation with RBF and Functions for Matrix builds
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
    # "Plotting points"
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

    
 # Basis Evaluation (which basis [i], where to evaluate [x, or a vector of x's], vector of centers for i to act on)
#  OLD WAY
def basis(x,xc):                     # RADIAL Basis i evaluated at x
    return ( exp(- beta * beta * (x-xc)*(x-xc) ) )
def basis_prime(x,xc):
    return( 2 * beta * beta * (xc-x) * exp( -beta * beta * (x-xc) * (x-xc) ) )
def basis_double_prime(x,xc):
    return( 2 * beta * beta * exp( -beta * beta * (x-xc)*(x-xc) ) * (2*beta*beta*(x-xc)*(x-xc) - 1) )
#def basis(x,xc):                    # RADIAL Basis i evaluated at x
#    beta_2 = beta * beta
#    xntr = (x-xc)
#    xntr_2 = xntr*xntr
#    return ( exp(- beta_2 * xntr_2 ) )
#def basis_prime(x,xc):              # Differential of Radial Basis evaluated at x
#    beta_2 = beta * beta
#    xntr = (x-xc)
#    xntr_2 = xntr*xntr
#    return( -2 * beta_2 * xntr * exp( - beta_2 * xntr_2 ) )
#def basis_double_prime(x,xc):       # Second Differential of Radial Basis evaluated at x
#    beta_2 = beta * beta
#    xntr = (x-xc)
#    xntr_2 = xntr*xntr
#    return( 2 * beta_2 * exp( -beta_2 * xntr_2 ) * (2*beta_2*xntr_2 - 1) )

    
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
        A[i,j]=p*basis_prime(xc[i],xc[j])-basis_double_prime(xc[i],xc[j])
# Build Matrix (Vector)  (adding one for initial/boundary condition)

g=np.zeros(np.size(xc)+2)                       # Same here with the "2"
#for i in range(np.size(xd)):
#    g[i] = f_double_prime(xd[i])

# Build out Matrix and vector for boundary equations
for j in range(np.size(xc)):                #  Boundary Condition 0
    A[np.size(xc)-1,j]=basis(xbound0,xc[j])
g[np.size(xc)-1]=f(xbound0)
for j in range(np.size(xc)):                #  Boundary Condition 1 (doesn't have to be at 1 and 0)
    A[np.size(xc),j]=basis(xbound1,xc[j])
g[np.size(xc)]=f(xbound1)

# Solve
a = np.linalg.lstsq(A,g)
coeff=a[0]      # Copy to coefficients

ans = f(xd)
#for i in range(np.size(xc)):
#    plot_rbf(i)    
plot_f()
plt.plot(xd,ans,'b')
