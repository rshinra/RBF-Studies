# -*- coding: utf-8 -*-
"""
Fokker Planck Attempt - 2d, using Dr, Meade's Template
"""
# Imports----------------------------------------------------------------------
import numpy as np
from math import exp
import matplotlib.pyplot as plt
# Functions--------------------------------------------------------------------
nxsample = 5
nysample = 5
nsample = nxsample * nysample
print "nsample: " + str(nsample)
nxeval = 5
nyeval = 5
neval = nxeval * nyeval
print "neval: " + str(neval)
D = .4
Xi = .2
Gam = .1
ranger = 2.0
print "ranger: " + str(ranger)
delstep = 2.0*ranger / (nxsample - 1)
print "delstep: " + str(delstep)
deleval = 2.0*ranger / (nxeval -1)
print "deleval: " + str(deleval)
eps = 1e-8
pdf = 1.0
xs = np.linspace(-ranger,ranger,nxsample)
ys = np.linspace(-ranger,ranger,nysample)
print xs,ys
X,Y =np.meshgrid(xs,ys,indexing='ij')
print X
print Y
def Phi(w,x1,x2,z1,z2):
    return np.exp(-w*w*((x1-z1)*(x1-z1)+(x2-z2)*(x2-z2)))
def bigA(w):
    return np.pi/(w*w)
def H(w,x1,x2,z1,z2):
    return ((w*w)*x2*(x1-z1))