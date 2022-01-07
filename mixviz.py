#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
import matplotlib.pyplot as plt
import h5py

# cartMap = True

def cartesianCellCenters(x,y,z):
    # Aliases to keep things short

    xc = 0.25*(x[:-1,:-1]+x[1:,:-1]+x[:-1,1:]+x[1:,1:])
    yc = 0.25*(y[:-1,:-1]+y[1:,:-1]+y[:-1,1:]+y[1:,1:])
    zc = 0.25*(z[:-1,:-1]+z[1:,:-1]+z[:-1,1:]+z[1:,1:])

    return xc,yc,zc

cmap = 'seismic'
cut = 0.5

hf = h5py.File('msphere.mix.h5','r')

x = np.array(hf['X'])

ntheta,nphi = x.shape

theta = np.linspace(0,np.pi,ntheta)
phi   = np.linspace(0,2*np.pi,nphi)

keys = np.array(list(hf.keys()))
nsteps = len(keys) - 2

print("nsteps=%d" %nsteps)

r = 1

th2D,phi2D = np.meshgrid(theta,phi,indexing='ij')

x = r * np.sin(th2D) * np.cos(phi2D)
y = r * np.sin(th2D) * np.sin(phi2D)
z = r * np.cos(th2D)

xc,yc,zc = cartesianCellCenters(x,y,z)

# Following are the keys in each step
# ['Average energy NORTH', 'Density NORTH', 'Field-aligned current NORTH', 'Hall conductance NORTH', 'IM Energy flux NORTH', 'IM Energy flux proton NORTH', 'IM average energy NORTH', 'IM average energy proton NORTH', 'Number flux NORTH', 'Pedersen conductance NORTH', 'Potential NORTH', 'Sound speed NORTH', 'Zhang average energy NORTH', 'Zhang number flux NORTH']

step = nsteps-100
plotKey = 'Field-aligned current'

datKey  = plotKey + ' NORTH'
dat = np.array(hf['Step#'+str(step)][datKey])
t   = hf['Step#' + str(step)].attrs.get('time',0.0)
datMax = (np.abs(dat)).max()

print(np.shape(dat),np.shape(th2D))


# fig = plt.figure(figsize=(12,12))

# ax = fig.add_subplot(projection='3d')

# sc = ax.scatter(xc,yc,zc,dat,c=dat,s=10*np.abs(dat),cmap=cmap,
#         vmin=-cut*datMax,vmax=cut*datMax)

# plt.colorbar(sc,fraction=0.06, pad=0.04,orientation='horizontal')
# ax.set_xlabel(r'$x$',fontsize=30)
# ax.set_ylabel(r'$y$',fontsize=30)
# ax.set_zlabel(r'$z$',fontsize=30)
# plt.tick_params(labelsize=20)

# if cartMap:
#     import cartopy.crs as ccrs
#     projection = ccrs.Mollweide()
#     ax  = fig.add_subplot(1,1,1,projection=projection)
#     phi = phi - np.pi
#     theta = np.pi/2 - theta
#     mesh = ax.pcolormesh(phi*180/np.pi,theta*180/np.pi,dat,cmap='seismic',transform=ccrs.PlateCarree(),vmin=-datMax,vmax=datMax,shading='auto')
#     plt.axis('off')
# else:
#     ax  = fig.add_subplot(1,1,1)
#     mesh = ax.pcolormesh(phi*180/np.pi,theta*180/np.pi,dat,cmap='seismic',vmin=-datMax,vmax=datMax,shading='auto')
#     ax.set_xlabel(r'$\phi$',fontsize=20)
#     ax.set_ylabel(r'$\theta$',fontsize=20)

# cax  = plt.colorbar(mesh,ax=ax)

# ax.set_title(plotKey + r'$, t=$%.3f' %t,fontsize=20)
# cax.ax.set_ylabel(plotKey,fontsize=15)

# if not cartMap:
#     plt.tight_layout()
# plt.show()

from mayavi import mlab

lut=eval('plt.cm.'+cmap+'(np.linspace(0,1,255))*255')

mlab.figure(size=(800, 800))
mesh_handle = mlab.mesh(xc,yc,zc, scalars=dat)
mesh_handle.module_manager.scalar_lut_manager.lut.table = lut
#    mesh_handle.module_manager.scalar_lut_manager.reverse_lut = True
mlab.draw()
mlab.show()