# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 09:47:30 2026

@author: Hermann
"""

import os
import numpy as np

dir0 = r"E:\Seg2Dat\Normandie\line10"
os.chdir(dir0)

times = np.loadtxt("picks.dat")
xs = (times[:,0]-1.)*4.
xr = times[:,1]-1.
d = np.abs(xs-xr)
t = times[:,2]
du = np.unique(d)
ndat = np.zeros_like(du, dtype=int)
tsum = np.zeros_like(du)
for i, tt in enumerate(t):
    n = np.where(d[i] == du)[0][0]
    ndat[n] += 1
    tsum[n] += tt
tav = tsum/ndat

with open("arrivals.dat", "w") as fo:
    for i, dd in enumerate(du):
        fo.write(f"{dd:2.0f} {tav[i]*1000:5.1f}\n")

