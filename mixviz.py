#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
import matplotlib.pyplot as plt
import h5py

cartMap = True

hf = h5py.File('msphere.mix.h5','r')

x = np.array(hf['X'])

ntheta,nphi = x.shape

theta = np.linspace(0,np.pi,ntheta)
phi   = np.linspace(0,2*np.pi,nphi)

keys = np.array(list(hf.keys()))
nsteps = len(keys) - 2

print("nsteps=%d" %nsteps)


# Following are the keys in each step
# ['Average energy NORTH', 'Density NORTH', 'Field-aligned current NORTH', 'Hall conductance NORTH', 'IM Energy flux NORTH', 'IM Energy flux proton NORTH', 'IM average energy NORTH', 'IM average energy proton NORTH', 'Number flux NORTH', 'Pedersen conductance NORTH', 'Potential NORTH', 'Sound speed NORTH', 'Zhang average energy NORTH', 'Zhang number flux NORTH']

step = nsteps-100
plotKey = 'Potential'

datKey  = plotKey + ' NORTH'
dat = np.array(hf['Step#'+str(step)][datKey])
t   = hf['Step#' + str(step)].attrs.get('time',0.0)
datMax = (np.abs(dat)).max()

fig = plt.figure(figsize=(10,10))

if cartMap:
    import cartopy.crs as ccrs
    projection = ccrs.Mollweide()
    ax  = fig.add_subplot(1,1,1,projection=projection)
    phi = phi - np.pi
    theta = np.pi/2 - theta
    mesh = ax.pcolormesh(phi*180/np.pi,theta*180/np.pi,dat,cmap='seismic',transform=ccrs.PlateCarree(),vmin=-datMax,vmax=datMax,shading='auto')
    plt.axis('off')
else:
    ax  = fig.add_subplot(1,1,1)
    mesh = ax.pcolormesh(phi*180/np.pi,theta*180/np.pi,dat,cmap='seismic',vmin=-datMax,vmax=datMax,shading='auto')
    ax.set_xlabel(r'$\phi$',fontsize=20)
    ax.set_ylabel(r'$\theta$',fontsize=20)

cax  = plt.colorbar(mesh,ax=ax)

ax.set_title(r'$t=$%.3f' %t,fontsize=20)
cax.ax.set_ylabel(plotKey,fontsize=15)

if not cartMap:
    plt.tight_layout()
plt.show()
