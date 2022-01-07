#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
import matplotlib.pyplot as plt
import h5py

def cartesianCellCenters(x,y,z):
    # Aliases to keep things short

    xc = 0.25*(x[:-1,:-1]+x[1:,:-1]+x[:-1,1:]+x[1:,1:])
    yc = 0.25*(y[:-1,:-1]+y[1:,:-1]+y[:-1,1:]+y[1:,1:])
    zc = 0.25*(z[:-1,:-1]+z[1:,:-1]+z[:-1,1:]+z[1:,1:])

    return xc,yc,zc

#cartMap = True

cmap='seismic'
cut= 0.5

hf = h5py.File('msphere.gam.h5','r')

x = np.array(hf['X'])
y = np.array(hf['Y'])
z = np.array(hf['Z'])

r     = np.sqrt(x**2 + y**2 + z**2)
theta = np.arccos(z/r)
phi   = np.arctan2(x,y)

keys = np.array(list(hf.keys()))
nsteps = len(keys) - 8

print("nsteps=%d" %nsteps)


# Following are the keys in each step
#['Bx', 'By', 'Bz', 'Cs', 'D', 'DivB', 'DivdB', 'Ex', 'Ey', 'Ez', 'Jx', 'Jy', 'Jz', 'P', 'Pb', 'Vx', 'Vy', 'Vz']

step = nsteps-100
datKey  = 'Vx'

Jx = np.array(hf['Step#'+str(step)]['Jx'])[...,0]
Jy = np.array(hf['Step#'+str(step)]['Jy'])[...,0]
Jz = np.array(hf['Step#'+str(step)]['Jz'])[...,0]

Bx = np.array(hf['Step#'+str(step)]['Bx'])[...,0]
By = np.array(hf['Step#'+str(step)]['By'])[...,0]
Bz = np.array(hf['Step#'+str(step)]['Bz'])[...,0]

x2d = x[...,0]
y2d = y[...,0]
z2d = z[...,0]

xc,yc,zc = cartesianCellCenters(x2d,y2d,z2d)

rc = np.sqrt(xc**2 + yc**2 + zc**2)

Jr = (Jx*xc + Jy*yc + Jz*zc)/rc
Br = (Bx*xc + By*yc + Bz*zc)/rc

dat = Br
t   = hf['Step#' + str(step)].attrs['time']
datMax = (np.abs(dat)).max()


# fig = plt.figure(figsize=(12,12))
# ax = fig.add_subplot(projection='3d')

# sc = ax.scatter(xc,yc,zc,dat,c=dat,s=10*np.abs(dat),cmap=cmap,
#         vmin=-cut*datMax,vmax=cut*datMax)

# plt.colorbar(sc,fraction=0.06, pad=0.04,orientation='horizontal')
# ax.set_xlabel(r'$x$',fontsize=30)
# ax.set_ylabel(r'$y$',fontsize=30)
# ax.set_zlabel(r'$z$',fontsize=30)
# plt.tick_params(labelsize=20)
# ax.set_title(r'FAC',fontsize=30)
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

# ax.set_title(r'$t=$%.3f' %t,fontsize=20)
# cax.ax.set_ylabel(datKey,fontsize=15)

# if not cartMap:
#     plt.tight_layout()

from mayavi import mlab

lut=eval('plt.cm.'+cmap+'(np.linspace(0,1,255))*255')

mlab.figure(size=(800, 800))
mesh_handle = mlab.mesh(xc,yc,zc, scalars=dat)
mesh_handle.module_manager.scalar_lut_manager.lut.table = lut
#    mesh_handle.module_manager.scalar_lut_manager.reverse_lut = True
mlab.draw()
mlab.show()