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
def u_hat(b,x,c_v,c_k):
    y = np.zeros(np.size(x))
    for i in range(np.size(x)):
        for j in range(np.size(c_v)):
            y[i] = y[i]+c_k[j]*basis(b,x[i],c_v[j])
    return y
def error(b,x,c_v,c_k):
    err = 0
    y = f(x)
    errv = np.abs(y-u_hat(b,x,c_v,c_k))
    err = np.sum(errv)    
    return err
 # Plot Final Function Approximation
def plot_f(ax,b,x,c_v,c_k):
    ax.plot(x,u_hat(b,x,c_v,c_k),'k')
    ax.plot(c_v,u_hat(b,c_v,c_v,c_k),'ko')
def plot_f_error(ax,b,x,c_v,c_k):
    y = f(x)
    ax.plot(x,y-u_hat(b,x,c_v,c_k),'k')    
# Plot Particular Basis Function
def plot_rbf(ax,i,b,x,c_v,c_k):
    y = np.zeros(np.size(x))
    for i in range(np.size(x)):
        y[i] = c_k[i]*basis(b,x[i],c_v[i])
    ax.plot(x,y,'r')
#------------------------------------------------------------------------------
# Main-------------------------------------------------------------------------
#Global variables--------------------------------------------------------------
  # Input parameters for RBF Level Ia
    # Basis center Locations / Number of Bases
xbound0 = -2.
xbound1 = 2.
xc=np.linspace(xbound0,xbound1,33)
    # Width parameter (How Fat your basis functions are)
betamin = 1.
betamax = 20.
betaspace = np.linspace(betamin,betamax,(betamax-betamin)*10+1)
    # "Evaluation points" (for function plotting)
xd = np.linspace(xbound0,xbound1,501)
  # Missing Parameter for RBF Level Ia
#coeff = np.ones(np.size(xc)) 
# Optimization Parameters
plt.clf
minerr=100
lam = np.zeros(np.size(betaspace))
SumC = np.zeros(np.size(betaspace))
Errr = np.zeros(np.size(betaspace))
for i in range(np.size(betaspace)):
    beta = betaspace[i]    
    H = Build_A(beta,xc)
    g = Build_b(np.size(xc))                       # Same here with the "2"
    # Boundary equations
    for j in range(np.size(xc)):                #  Boundary Condition 0, integral, i.e. all ck must add up to that
        H[np.size(xc),j]=np.sqrt(np.pi)/beta
    g[np.size(xc)]=1
    
    # Solve
    a = np.linalg.lstsq(H,g)
    ck=a[0]      # Copy to coefficients
    if (np.abs(ck[np.size(xc)]) > .000001):
      lam[i] = ck[np.size(xc)]
    else:
      lam[i] = 0
    SumC[i] = sum(H[np.size(xc)]*ck)
    Errr[i] = error(beta,xc,xc,ck)    
    if(Errr[i] < minerr):
        betafinal = beta
        ckfinal = ck
        lamfinal = lam[i]
        SCfinal = SumC[i]
        minerr = error(beta,xc,xc,ck)
ans=f(xd)    
print(betafinal)
print(lamfinal)
print(SCfinalg)
ax1 = plt.subplot(221)
ax2 = plt.subplot(223)
ax3 = plt.subplot(122)
plot_f_error(ax3,betafinal,xc,xc,ckfinal)
plot_f(ax1,betafinal,xd,xc,ckfinal)
ax1.plot(xd,ans,'b')
ax2.plot(betaspace,Errr,'b')
