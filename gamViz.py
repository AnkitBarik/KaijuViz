#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from kaipy.gamera.gampp import *
from kaipy.gamera.magsphere import *
from kaipy.remix.remix import *
import kaipy.kaiH5 as kaiH5
from astropy.time import Time
from scipy.interpolate import griddata 

def get2Dgrid(s,phi,Jr_sur_NH):
    
    s2D = np.zeros_like(Jr_sur_NH[...,0])
    p2D = np.zeros_like(Jr_sur_NH[...,0])

    sInterp = np.linspace(0.,s.max(),s2D.shape[0])
    phiInterp = np.linspace(0.,2*np.pi,s2D.shape[1])
    
    for k,sk in enumerate(sInterp):
        s2D[k,:] = sk
    for k,pk in enumerate(phiInterp):
        p2D[:,k] = pk
    
    x2D = s2D * np.cos(p2D)
    y2D = s2D * np.sin(p2D)

    return x2D, y2D

def interp_dat(xx,yy,zz,nk,hem,dat):
    
    half = int(nk/2)
    if hem.upper()=="NORTH":
        xx_1 = xx[...,:half]
        yy_1 = yy[...,:half]
        zz_1 = zz[...,:half]
        dat  =dat[...,:half]

    elif hem.upper() == "SOUTH":
        xx_1 = xx[...,half:]
        yy_1 = yy[...,half:]
        zz_1 = zz[...,half:]
        dat  =dat[...,half:]

    r = np.sqrt(xx_n**2 + yy_n**2 + zz_n**2)

    mask = r <= r[0,0,0]

    x_mask = xx_1[mask]
    y_mask = yy_1[mask]

    s = np.sqrt(x_mask**2 + y_mask**2)
    phi = np.arctan2(y_mask,x_mask)

    xx2,yy2 = get2Dgrid(s,phi,Jr_sur_NH)
#theta = np.arctan2(s,zz_n)

    dat_new = griddata((x_mask,y_mask), Jr_sur_NH[mask], (xx2,yy2),method='linear')




tstep = 120

gIn = GamsphPipe(fdir='.',ftag='msphere')

Jx = gIn.GetVar("Jx",tstep)
Jy = gIn.GetVar("Jy",tstep)
Jz = gIn.GetVar("Jz",tstep)

#Bx = gIn.GetVar("Bx",tstep)
#By = gIn.GetVar("By",tstep)
#Bz = gIn.GetVar("Bz",tstep)

xx = gIn.X[:-1,:-1,:-1]
yy = gIn.Y[:-1,:-1,:-1]
zz = gIn.Z[:-1,:-1,:-1]

rnorm = np.sqrt(xx**2 + yy**2 + zz**2)
irnorm = 1./rnorm

Jr_surface = irnorm * (Jx*xx + Jy*yy + Jz*zz)

xScl = 1000e3 #m

bScl = 4.58e-9 #nT

mu0 = 4*np.pi * 1e-7

JScl = bScl/xScl /mu0 * 1e6 #micro A/m^2

dat_new *= JScl

plt.figure(figsize=(10.2,10))
#plt.axis('equal')
plt.contourf(xx2,yy2,dat_new,50,cmap='RdBu_r')
plt.colorbar()

datMax = max(np.nanmax(dat_new),np.nanmin(dat_new))
plt.clim(-datMax,datMax)

plt.axis('equal')
plt.axis('off')
plt.tight_layout()
#plt.show()
plt.savefig('gamviz.png',dpi=150)

#plt.figure(figsize=(16,9))




#plt.contourf(xx[:,:,0],yy[:,:,0],Jr_surface[:,:,0],50,cmap='RdBu_r')
#plt.contourf(xx[:,:,-1],yy[:,:,-1],Jr_surface[:,:,-1],50,cmap='RdBu_r')
#
#plt.axis('equal')
#plt.axis('off')
#plt.show()            
