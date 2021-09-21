#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
import matplotlib.pyplot as plt
import h5py

cartMap = False

hf = h5py.File('msphere.gam.h5','r')

x = np.array(hf['X'])
y = np.array(hf['Y'])
z = np.array(hf['Z'])

r     = np.sqrt(x**2 + y**2 + z**2)
theta = np.arccos(z/r)[...,0]
phi   = np.arctan2(x,y)[...,0]

keys = np.array(list(hf.keys()))
nsteps = len(keys) - 8

print("nsteps=%d" %nsteps)


# Following are the keys in each step
#['Bx', 'By', 'Bz', 'Cs', 'D', 'DivB', 'DivdB', 'Ex', 'Ey', 'Ez', 'Jx', 'Jy', 'Jz', 'P', 'Pb', 'Vx', 'Vy', 'Vz']

step = nsteps-100
datKey  = 'Vx'


dat = np.array(hf['Step#'+str(step)][datKey])[...,0]
t   = hf['Step#' + str(step)].attrs['time']
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
cax.ax.set_ylabel(datKey,fontsize=15)

if not cartMap:
    plt.tight_layout()
plt.show()
