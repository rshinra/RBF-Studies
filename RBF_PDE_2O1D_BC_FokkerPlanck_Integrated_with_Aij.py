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
# Functions--------------------------------------------------------------------
# Function for evaluation of error:  f(x) -> 
def f(x):
   return .455 * np.exp(-10*x*x*x*x + 5*x*x)
 # Function Derivative for "test Data" f'(x) -> 
def f_prime(x):
    return .455 * np.exp(-10*x*x*x*x + 5*x*x) * (-40*x*x*x + 10*x)
 # Function Second Derivative for "test Data" f''(x) -> 1
def f_double_prime(x):
    return .455 * np.exp(-10*x*x*x*x + 5*x*x) * (-120*x*x + 10 + (-40*x*x*x + 10*x)*(-40*x*x*x + 10*x))

# GAMMA function (r)
def gamma(x):
    return -40*x*x*x + 10*x
def gamma_prime(x):
    return -120*x*x + 10
    
 # Basis Evaluation
def basis(b,x,c):                     # RADIAL Basis i evaluated at x
    b_2 = b*b
    x_c = x-c
    x_c_2 = x_c*x_c
    return ( exp(-b_2 * x_c_2 ) )    
 # Differential of basis; Recall, A Differential of a SUM is a SUM of Differentials
def basis_prime(b,x,c):
    b_2 = b*b
    x_c = x-c
    x_c_2 = x_c*x_c
    return( -2 * b_2 * x_c * exp( -b_2 * x_c_2 ) )
 # Second differential of basis; Recall, A Differential of a SUM is a SUM of Differentials
def basis_double_prime(b,x,c):
    b_2 = b*b
    x_c = x-c
    x_c_2 = x_c*x_c
    return( 2*b_2*exp( -b_2 * x_c_2 ) * (2*b_2*x_c_2 - 1) )
# Matrixes for the problem
def Build_A(b,c_v):               # Matrix H (Ax=b, Hc=g)
    A = np.zeros((np.size(c_v)+1,np.size(c_v)+1))       # Why is there a 1 here?  (Integral Boundary Condition)
    for j in range(np.size(c_v)):     # J is columns of A
        for i in range(np.size(c_v)):    # I is rows of A
            A[i,j]=Aij(b,c_v[i],c_v[j])
    for i in range(np.size(c_v)):
        A[i,np.size(c_v)]=1.         # for Lagrange Multiplier
    return A
def Build_b(l):
    b = np.zeros(np.size(xc)+1)
    return b
def Aij(b,x,c):
    return basis_double_prime(beta,x,c)-basis_prime(beta,x,c)*gamma(x)-basis(beta,x,c)*gamma_prime(x)
def u_hat(b,x,c_v):
    y = np.zeros(np.size(x))
    for i in range(np.size(x)):
        for j in range(np.size(c_v)):
            y[i] = y[i]+coeff[j]*basis(b,x[i],c_v[j])
    return y
    
 # Plot Final Function Approximation
def plot_f(b,x,c_v):
    plt.plot(x,u_hat(b,x,c_v),'k')
    plt.plot(c_v,u_hat(b,c_v,c_v),'ko')
    
# Plot Particular Basis Function
def plot_rbf(i,b,x,c_v,c_k):
    y = np.zeros(np.size(x))
    for i in range(np.size(x)):
        y[i] = c_k[i]*basis(b,x[i],c_v[i])
    plt.plot(xd,y,'r')
#------------------------------------------------------------------------------
# Main-------------------------------------------------------------------------
#Global variables--------------------------------------------------------------
  # Input parameters for RBF Level Ia
    # Basis center Locations / Number of Bases
xbound0 = -2
xbound1 = 2
xc=np.linspace(xbound0,xbound1,33)
    # Width parameter (How Fat your basis functions are)
beta = 3
    # "Evaluation points" (for function plotting)
xd = np.linspace(xbound0,xbound1,501)
  # Missing Parameter for RBF Level Ia
#coeff = np.ones(np.size(xc)) 
plt.clf()   #clear plots
H = Build_A(beta,xc)
g = Build_b(np.size(xc))                       # Same here with the "2"
# Boundary equations
for j in range(np.size(xc)):                #  Boundary Condition 0, integral, i.e. all ck must add up to that
    H[np.size(xc),j]=np.sqrt(np.pi)/beta
g[np.size(xc)]=1

# Solve
a = np.linalg.lstsq(H,g)
coeff=a[0]      # Copy to coefficients

ans = f(xd)
#for i in range(np.size(xc)):
#    plot_rbf(i)    
plot_f(beta,xd,xc)
plt.plot(xd,ans,'b')
#print coeff[np.size(xc)]
#print sum(A[201]*coeff)
