#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
import sys
from kaipy.gamera.gampp import *
from kaipy.remix.remix import *
import kaipy.kaiH5 as kaiH5
from astropy.time import Time
import cartopy.crs as ccrs

def get_latlong(rem):

    x = rem.ion['X']
    y = rem.ion['Y'] 
    phi_nh = rem.ion['THETA'][1:,:]

    phi_nh -= phi_nh.min()

    s = np.sqrt(x**2 + y**2)
    s /= s.max()

    theta_nh = np.arcsin(s)
    theta_nh = theta_nh[1:,:]

    theta_sh = theta_nh + np.pi/2
    phi_sh = phi_nh[:,::-1]
    
    theta =  np.pi/2 - np.vstack([theta_nh,theta_sh])
    phi   = np.vstack([phi_nh,phi_sh]) - np.pi

    return theta, phi

def symmetrize(dat1):
    
    dat = np.zeros([dat1.shape[0],dat1.shape[1]+1])
    dat[:,:-1] = dat1
    dat[:,-1] = dat1[:,0]

    return dat

    
def get_dat(rem,datType):
    
    if datType.lower()=='fac':

        fac_nh1 = rem.ion['Field-aligned current NORTH']
        fac_nh = symmetrize(fac_nh1)
        

        fac_sh1 = rem.ion['Field-aligned current SOUTH']
        fac_sh = symmetrize(fac_sh1)
        
        fac   = np.vstack([fac_nh,fac_sh])

        return fac

    elif datType.lower()=='potential':

        pot_nh1 = rem.ion['Potential NORTH']
        pot_nh = symmetrize(pot_nh1)
        pot_sh1 = rem.ion['Potential SOUTH']
        pot_sh = symmetrize(pot_sh1)

        pot = np.vstack([pot_nh,pot_sh])

        return pot
        

remixFile = 'msphere.mix.h5'

## Get time stamps

nsteps,sIds=kaiH5.cntSteps(remixFile)
T=kaiH5.getTs(remixFile,sIds,aID='MJD')

print("nsteps=",nsteps)

t = Time(T, format='mjd').iso

rem = remix(remixFile,1)
theta, phi = get_latlong(rem)

plotcrs = ccrs.Mollweide()


for tstep in range(1,2):

    rem = remix(remixFile,tstep)
    timeTitl = str.split(t[tstep]," ")[1]

    fac = get_dat(rem,'fac')
    pot = get_dat(rem,'potential')

# Plotting remix data
    fig = plt.figure(figsize=(10,12))

    ax1 = plt.subplot(2,1,1,projection=plotcrs)
#    divnorm = colors.TwoSlopeNorm(vmin=fac.min(), vcenter=0, vmax=fac.max())
    cont = ax1.contourf(phi*180/np.pi,theta*180/np.pi,fac,100,transform=ccrs.PlateCarree(),cmap='RdBu_r')#,norm=divnorm)

    plt.colorbar(cont,ax=ax1)
    ax1.set_title("FAC",fontsize=30)

    ax2 = plt.subplot(2,1,2,projection=plotcrs)
#    divnorm = colors.TwoSlopeNorm(vmin=pot.min(), vcenter=0, vmax=pot.max())
    cont = ax2.contourf(phi*180/np.pi,theta*180/np.pi,pot,100,transform=ccrs.PlateCarree(),cmap='PRGn')#,norm=divnorm)

    plt.colorbar(cont,ax=ax2)
    ax2.set_title("Potential",fontsize=30)

    plt.gcf().suptitle("t = %s" %timeTitl,fontsize=20)

    plt.tight_layout()
    plt.show()
#    plt.savefig('mixpic/img%03d.png' %tstep, dpi=150)
#    print("%d/%d" %(tstep,nsteps))
#    plt.close()

