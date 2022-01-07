#!/usr/bin/env python3

from pylab import *
from libaod import *
from magic.libmagic import *

r_surface=2.4e6 # metres

eps = 0.

r_mer = 2.4

r = 1e6*chebgrid(65,r_mer,5*r_mer) #metres

#r = 1e6*linspace(2.4,4.8,65) # metres

theta = linspace(0+eps,np.pi-eps,128)
phi = linspace(0,2*np.pi,256)

eta = r.min()/r.max()

r3D, th3D, p3D, x, y, z, s = get_grid(r,theta,phi,len(r),len(theta),len(phi))

BrSPH = getB_SPH_mer(r3D,p3D,th3D,r_surface)

mag_moment_dip  = -190e-9 # Tesla
mag_moment_quad = -74.6e-9 # Tesla

m_dip  = r_surface**3 * mag_moment_dip
m_quad = r_surface**4 * mag_moment_quad

Br1 = BCart(x,y,z,m_dip, m_quad, th3D,p3D)
Br2 = BVecPot(x,y,z,m_dip,m_quad,r3D,th3D,p3D,eta,theta)

err = abs(Br1[0,...] - BrSPH[0,...]) + abs(Br2[0,...] - BrSPH[0,...])

Valf = getValf(x,y,z,m_dip, m_quad, th3D, p3D)

cmap = 'seismic'

radLevel = 0

x2d = x[...,radLevel]
y2d = y[...,radLevel]
z2d = z[...,radLevel]

dat = Br1[...,radLevel]

from mayavi import mlab

lut=eval('plt.cm.'+cmap+'(np.linspace(0,1,255))*255')

mlab.figure(size=(800, 800))
mesh_handle = mlab.mesh(x2d,y2d,z2d, scalars=dat)
mesh_handle.module_manager.scalar_lut_manager.lut.table = lut
#    mesh_handle.module_manager.scalar_lut_manager.reverse_lut = True
mlab.draw()
mlab.show()

# nplots = 5

# figure(figsize=(nplots*5.1,10))

# subplot(1,nplots,1)
# merContour(r,theta,BrSPH[30,...]*1e9)
# title("Spherical Harmonics")

# subplot(1,nplots,2)
# merContour(r,theta,Br1[30,...]*1e9)
# title("Magnetic moments")

# subplot(1,nplots,3)
# merContour(r,theta,Br2[30,...]*1e9)
# title("Vector Potential")

# errPlotIdx = 2

# subplot(1,nplots,4)
# merContour(r,theta[errPlotIdx:-errPlotIdx],err[errPlotIdx:-errPlotIdx,:]*1e9)
# title("Sum of errors")

# subplot(1,nplots,5)
# merContour(r,theta,Valf[0,...])
# title(r"Alfvén speed")

# tight_layout()
# show()
#savefig('BrBenchmark.pdf',dpi=400)