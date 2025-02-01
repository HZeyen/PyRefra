# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 09:38:35 2023

@author: Hermann Zeyen
         University Paris-Saclay

Program reads file "velocities.dat" written by program PyRefra after tomographic
inversion.

User should set the following parameters:
    dir0 : The folder where file velocities.dat is located
    dx, dy : grid size for interpolation.
    v_search : the velocity to be extracted in m/s
    title : Title for plot (the chosen velocity is added automatically to this text)

"""

import numpy as np
import os
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Set user inputs
dir0 = r"E:/Seg2Dat/Barbeau/20230628-19.46"
dx = 1.
dy = 0.2
v_search = 1500.
title = "Barbeau"

# Read data file. Sign of Y coordinates is inverted since PyGimli uses positive
#      upward coordinate system
os.chdir(dir0)
data = np.loadtxt("velocities.dat")
x = data[:,0]
y = -data[:,1]
v = data[:,2]

# Extract grid limits as multiples of dx and dy. Adding/subtracting dx or dy
#    is eventually necessary to avoid extrapolation
xmn = np.round(x.min()/dx,0)*dx
if xmn < x.min():
    xmn += dx
xmx = np.round(x.max()/dx,0)*dx
if xmx > x.max():
    xmx -= dx
ymn = np.round(y.min()/dy,0)*dy
if ymn < y.min():
    ymn += dy
ymx = np.round(y.max()/dy,0)*dy
if ymx > y.max():
    ymx -= dy

# Interpolate velocities to regular grid
X_unique = np.arange(xmn,xmx+dx/2,dx)
Y_unique = np.arange(ymn,ymx+dy/2,dy)
X, Y = np.meshgrid(X_unique, Y_unique)  # 2D grid for interpolation
interp = LinearNDInterpolator(list(zip(x, y)), v)
V = interp(X, Y)

# Go through every column and find first occurrence of v > desired velocity
# Then calculate the precise position through linear interpolation
z_iso = []
x_iso = []
for i in range(len(X_unique)):
    n = np.where(V[:,i] > v_search)[0]
    if len(n) > 0:
        j = n[0]
        x_iso.append(X_unique[i])
        z_iso.append((Y_unique[j]-Y_unique[j-1])/(V[j,i]-V[j-1,i])*(v_search-V[j-1,i])+\
                      Y_unique[j-1])
    else:
        continue

# Plot velocities and velocity isoline
fig,ax = plt.subplots(1,1, figsize=(12,8))
pl = ax.imshow(V, cmap="rainbow",extent=[xmn,xmx,ymx,ymn])
ax.plot(x_iso,z_iso,"k")
ax.set_xlabel("Distance[m]")
ax.set_ylabel("Depth [m]")
ax.set_title(f"{title} Isoline v={v_search:0.0f} m/s")
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.2)
cb = fig.colorbar(pl, cax=cax, label="velocity [m/s]")
fig.savefig(f"isoline_{v_search:0.0f}.png")

# Save isoline coordinates to file
with open(f"isoline_{v_search:0.0f}.txt","w") as fo:
    fo.write(f"x, z for isoline v={v_search:0.0f} m/s\n")
    for i in range(len(x_iso)):
        fo.write(f"{x_iso[i]:0.1f}  {z_iso[i]:0.1f}\n")

