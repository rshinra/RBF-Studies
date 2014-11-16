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
nxsample = 33
nysample = 1
nsample = nxsample * nysample
nxeval = 33
nyeval = 1
neval = nxeval*nyeval
dee = .1
gam = .5
mu = 0.0
plusminus = 2.0
delstep = 2.0 * plusminus / (nxsample - 1)
eps = 1e-8
pdf = 1
xsvec = np.linspace(-plusminus,plusminus,nxsample)
print xsvec
evalcoords = xsvec
def Phi(w,x1,z1):
    return np.exp(-((x1-z1)*(x1-z1))*w*w)
def BigA(w):
    return np.sqrt(np.pi/(w*w))
etaa1=40.0
etaa2=10.0
def H(w,x1,z1):
    return (-2.0*(w*w) + 4.0*(w*w*w*w)*(x1-z1)*(x1-z1) + (3.0*40.0*(x1*x1)-10) - 2.0*(w*w)*(x1-z1)*(40*x1*x1*x1-10.0*x1))*np.exp(-((x1-z1)*(x1-z1))*w*w)
print H(1.76,evalcoords[0],xsvec[0])
def kahalu(w,i,j):
    return H(w,evalcoords[i], xsvec[j])
rhsvector = np.zeros(nsample+1)
rhsvector[nsample] = pdf
print rhsvector
oldthesigma = 5.1
thesigma = 3.0
mauka = np.ones(nsample+1,nsample+1)
for i in range(nsample):
    for j in range(nsample):
        mauka[i,j] = kahalu[thesigma,i,j]
for i in range nsam
            