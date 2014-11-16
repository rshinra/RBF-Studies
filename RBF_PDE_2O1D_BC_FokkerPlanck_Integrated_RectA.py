# -*- coding: utf-8 -*-
"""
Fokker Planck Attempt - WORKS, but with wrong size of matrix A (or H)
(Solving Differential Equations using Radial Basis functions)
Function to model:      r = -40x^3 + 10x | u = .455 e^[-10x^4+5x^2]
Differential equation:  (ru)' = u'' | Integral of u from -oo to oo is 1
David Norman
"""
# Imports----------------------------------------------------------------------
import numpy as np
from math import exp
import matplotlib.pyplot as plt
#Global variables--------------------------------------------------------------
  # Input parameters for RBF Level Ia
    # Basis center Locations / Number of Bases
xbound0 = -2
xbound1 = 2
xc=np.linspace(xbound0,xbound1,50)
    # Width parameter (How Fat your basis functions are)
#beta = np.linspace(1,20,39)
beta = 4
beta2 = beta*beta
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
A = np.zeros((np.size(xd)+1,np.size(xc)+1))       # Why is there a 1 here?  (number of bc's)
for j in range(np.size(xc)):     # J is columns of A
    for i in range(np.size(xd)):    # I is rows of A
        A[i,j]=basis_prime(j,xd[i])-basis(j,xd[i])*gamma(xd[i])
for i in range(np.size(xd)):
    A[i,50]=1.
# Build Matrix (Vector)  (adding one for initial/boundary condition)

g=np.zeros(np.size(xd)+1)                       # Same here with the "2"
#for i in range(np.size(xd)):
#    g[i] = f_double_prime(xd[i])

# Build out Matrix and vector for boundary equations
for j in range(np.size(xc)):                #  Boundary Condition 0
    A[np.size(xd),j]=np.sqrt(np.pi)/beta
g[np.size(xd)]=1

# Solve
a = np.linalg.lstsq(A,g)
coeff=a[0]      # Copy to coefficients

ans = f(xd)
#for i in range(np.size(xc)):
#    plot_rbf(i)    
plot_f()
plt.plot(xd,ans,'b')
#print coeff[np.size(xc)]
#print sum(A[201]*coeff)
