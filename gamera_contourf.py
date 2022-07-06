#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import h5py

def cartesianCellCenters(x,y,z):

    xc = 0.25*(x[:-1,:-1]+x[1:,:-1]+x[:-1,1:]+x[1:,1:])
    yc = 0.25*(y[:-1,:-1]+y[1:,:-1]+y[:-1,1:]+y[1:,1:])
    zc = 0.25*(z[:-1,:-1]+z[1:,:-1]+z[:-1,1:]+z[1:,1:])

    return xc,yc,zc

def tileAndFix(dat):
    '''
    Tile around for phi=2\pi and append it
    '''
    dat_r = np.zeros([dat.shape[0]+1,dat.shape[1]])
    dat_r[:-1,:] = dat
    dat_r[-1,:] = dat[0,:]

    return dat_r

def noonMidnightSlice(x,y,z,dat,nphi,cmap='seismic'):

    '''
    Noon-midnight meridian
    '''

    idx1 = nphi//4 #pi/2
    idx2 = 3*nphi//4 #3\pi/2

    fig,ax = plt.subplots(figsize=(10,10))

    x2d  = x[idx1,...]
    y2d  = y[idx1,...]
    z2d  = z[idx1,...]
    dat1 = dat[idx1,...]

    xc1,yc1,zc1 = cartesianCellCenters(x2d,y2d,z2d)

    ax.contourf(xc1,zc1,dat1,100,cmap=cmap)

    x2d  = x[idx2,...]
    y2d  = y[idx2,...]
    z2d  = z[idx2,...]
    dat2 = dat[idx2,...]

    xc2,yc2,zc2 = cartesianCellCenters(x2d,y2d,z2d)

    ax.contourf(xc2,zc2,dat2,100,cmap=cmap)

    ax.invert_xaxis()
    ax.set_aspect('equal')
    ax.axis('off')

    return ax

def dawnDuskSlice(x,y,z,dat,ntheta,cmap='seismic'):

    '''
    Dawn Dusk slice
    '''

    idx = ntheta//2

    fig,ax = plt.subplots(figsize=(10,10))

    x2d     = x[:,idx,:]
    y2d     = y[:,idx,:]
    z2d     = z[:,idx,:]
    dat2d   = dat[:,idx,:]

    xc,yc,zc = cartesianCellCenters(x2d,y2d,z2d)

    yc_r  = tileAndFix(yc)
    zc_r  = tileAndFix(zc)
    dat_r = tileAndFix(dat2d)

    ax.contourf(yc_r,zc_r,dat_r,100,cmap=cmap)
    ax.set_aspect('equal')
    ax.axis('off')

hf = h5py.File('msphere.gam.h5','r')

x = np.array(hf['X'])
y = np.array(hf['Y'])
z = np.array(hf['Z'])

r = np.sqrt(x**2 + y**2 + z**2)
theta = np.arccos(x/r)
phi = np.arctan2(y,z)

(nphi,ntheta,nr) = x.shape

keys = np.array(list(hf.keys()))
nsteps = len(keys) - 8

print("nsteps=%d" %nsteps)


# Following are the keys in each step
#['Bx', 'By', 'Bz', 'Cs', 'D', 'DivB', 'DivdB', 'Ex', 'Ey', 'Ez', 'Jx', 'Jy', 'Jz', 'P', 'Pb', 'Vx', 'Vy', 'Vz']

step = nsteps-10
datKey  = 'Cs'
t   = hf['Step#' + str(step)].attrs['time']
dat = np.array(hf['Step#'+str(step)+'/'+datKey])
cmap = 'afmhot'

noonMidnightSlice(x,y,z,dat,nphi,cmap=cmap)
dawnDuskSlice(x,y,z,dat,ntheta,cmap=cmap)

plt.show()