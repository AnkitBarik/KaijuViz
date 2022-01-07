#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import h5py
import cartopy.crs as ccrs

plt.style.use('dark_background')

def getCellCenters(x,y,z):
    xc = (1./8.) *   ( x[:-1,:-1,:-1]
                    +  x[1:,:-1,:-1]  + x[1:,1:,:-1]
                    +  x[:-1,1:,:-1]  + x[:-1,1:,1:]
                    +  x[:-1,:-1,1:]  + x[1:,:-1,1:]
                    +  x[1:,1:,1:] )
    yc = (1./8.) *   ( y[:-1,:-1,:-1]
                    +  y[1:,:-1,:-1]  + y[1:,1:,:-1]
                    +  y[:-1,1:,:-1]  + y[:-1,1:,1:]
                    +  y[:-1,:-1,1:]  + y[1:,:-1,1:]
                    +  y[1:,1:,1:] )
    zc = (1./8.) *   ( z[:-1,:-1,:-1]
                    +  z[1:,:-1,:-1]  + z[1:,1:,:-1]
                    +  z[:-1,1:,:-1]  + z[:-1,1:,1:]
                    +  z[:-1,:-1,1:]  + z[1:,:-1,1:]
                    +  z[1:,1:,1:] )

    return xc,yc,zc

def plotSurf(p2D,th2D,data,ax,levels=60,cmap='RdBu_r'):

    bmax = np.abs(data).max()

    divnorm = colors.TwoSlopeNorm(vmin=-bmax, vcenter=0, vmax=bmax)
    cs = np.linspace(-bmax,bmax,levels)


    lon2D = p2D #- np.pi (phi is already from -pi tp pi)
    lat2D = np.pi/2 - th2D

    print(np.abs(th2D).max())
    print(np.abs(p2D).max())

    cont = ax.contourf(lon2D*180/np.pi,lat2D*180/np.pi,data,cs,
        transform=ccrs.PlateCarree(),
        cmap=cmap,norm=divnorm,extend='both')

    cbar = plt.colorbar(cont,orientation='vertical',fraction=0.06, pad=0.04,ticks=[-bmax,0,bmax])

    ax.axis('equal')
    ax.axis('off')

    return ax, cbar

hf = h5py.File('msphere.volt.h5','r')

x = np.array(hf['X'])
y = np.array(hf['Y'])
z = np.array(hf['Z'])

x,y,z = getCellCenters(x,y,z)

r      = np.sqrt(x**2 + y**2 + z**2)
theta  = np.arccos(z/r)
phi    = np.arctan2(y,x)

keys = np.array(list(hf.keys()))

nsteps = len(keys) - 3

proj='Mollweide'
projection = eval('ccrs.'+proj+'()')
cmap = 'RdBu_r'

# Choose a data field from ['Cs', 'D', 'Ex', 'Ey', 'Ez', 'Jx', 'Jy', 'Jz', 'Vx', 'Vy', 'Vz', 'dBx', 'dBy', 'dBz', 'psi']

field = 'psi'
removeShell = False

cut = 0.5
step = nsteps - 100

dat = np.array(hf['Step#'+str(step)+'/'+field])

radLevel = 1

x   = x[...,radLevel]
y   = y[...,radLevel]
z   = z[...,radLevel]
dat = dat[...,radLevel]

# for step in [nsteps-1]:
#     dat = np.array(hf['Step#'+str(step)+'/'+field])

#     datMax = (np.abs(dat)).max()

#     fig = plt.figure(figsize=(12,12))
#     ax = fig.add_subplot(projection='3d')
#     if removeShell:
#         x   = x[...,:-1]
#         y   = y[...,:-1]
#         z   = z[...,:-1]
#         dat = dat[...,:-1]
#     sc = ax.scatter(x,y,z,dat,c=dat,s=1e2*np.abs(dat),cmap=cmap,vmin=-cut*datMax,vmax=cut*datMax)
#     plt.colorbar(sc,fraction=0.06, pad=0.04,orientation='horizontal')
#     ax.set_xlabel(r'$x$',fontsize=30)
#     ax.set_ylabel(r'$y$',fontsize=30)
#     ax.set_zlabel(r'$z$',fontsize=30)
#     plt.tick_params(labelsize=20)
#     ax.set_title(r'$\psi$',fontsize=30)

# plt.tight_layout()
# plt.show()

from mayavi import mlab

lut=eval('plt.cm.'+cmap+'(np.linspace(0,1,255))*255')

mlab.figure(size=(800, 800))
mesh_handle = mlab.mesh(x,y,z, scalars=dat)
mesh_handle.module_manager.scalar_lut_manager.lut.table = lut
#    mesh_handle.module_manager.scalar_lut_manager.reverse_lut = True
mlab.title('psi')
mlab.draw()
mlab.show()