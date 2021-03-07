#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from kaipy.gamera.gampp import *
from kaipy.remix.remix import *
import kaipy.kaiH5 as kaiH5

remixFile = 'msphere.mix.h5'

nsteps,sIds=kaiH5.cntSteps(remixFile)

maxdata = np.zeros(nsteps)
mindata = np.zeros(nsteps)

dataType = 'current'
#dataType = 'potential'

for tstep in range(1,nsteps):
    rem = remix(remixFile,tstep)

    rem.init_vars("NORTH")
    FAC_mix_NH = rem.variables[dataType]['data']

    rem.init_vars("SOUTH")
    FAC_mix_SH = rem.variables[dataType]['data']
    
    max_NH = FAC_mix_NH.max()
    min_NH = FAC_mix_NH.min()
    max_SH = FAC_mix_SH.max()
    min_SH = FAC_mix_SH.min()

    maxdata[tstep] = max(max_NH,max_SH)
    mindata[tstep] = min(min_NH,min_SH)


maxDat = maxdata.max()
minDat = mindata.min()

np.savetxt("mix_%s_lim.dat" %dataType,[minDat, maxDat])
