# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 10:00:46 2024

@author: Hermann Zeyen
         University Paris-Saclay

Program transforms a velocity-depth file obtained from tomographic inversion
with program PyRefra into RMS velocities vs TWT used with VISTA

The input file is usually calles velocities.dat. If not, the name has to be
modified in line "data = np.loadtxt..."
User should first set manually the working directory in line "dir0 = ..."

The program interpolates first the velocities onto a regular grid, calculates
then the zero-offset TWT and then the RMS velocities.

Results are plotted in three graphs and written into file velocities.vel
This file should be copied into the folder "MISC" of the VISTA project and
imported from within the 2D Velocity Analysis window.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata

sys_path = r"E:/Sources_2010/Python_programs"
if sys_path not in sys.path:
    sys.path.append(sys_path)
from Dialog_HZ import dialog

dir0 = r"E:/Seg2Dat/Fontaines-Salees/2023/Line1/20230429-23.21"
v_file = "velocities.dat"
x_cmp1 = -0.5
dx_cmp = 0.5
dy = 0.1
n_cmp_start = 30
n_cmp_end = 330
d_n_cmp = 30
n_dy = 10
tmax = 100
dx = dx_cmp

results, button = dialog(["Velocity file",
                          "X0 of CMP 1 [m]",
                          "Distance between CMPs [m]",
                          "Vertical grid step [m]",
                          "maximum time [ms]",
                          "\nFollowing CMP #: natural counting",
                          "First CMP to be treated",
                          "Last CMP to be treated",
                          "CMP step",
                          "Output every N vertical grid points"],\
                         ["e","e","e","e","e","l","e","e","e","e"],\
                         [v_file,x_cmp1,dx_cmp,dy,tmax,None,n_cmp_start,n_cmp_end,d_n_cmp,n_dy],\
                         "CMP information")
v_file = results[0]
x_cmp1 = float(results[1])
dx_cmp = float(results[2])
dy = float(results[3])
tmax = float(results[4])/1000.
n_cmp_start = int(results[6])
n_cmp_end = int(results[7])
d_n_cmp = int(results[8])
n_dy = int(results[9])

os.chdir(dir0)
data = np.loadtxt(v_file, dtype=float)
if len(v_file) == 0:
    print("Velocity file has no data, program stops.")
    sys.exit()
x = data[:,0]    
y = data[:,1]    
vel = data[:,2]
xmin = np.ceil(x.min()/dx)*dx
xmax = np.floor(x.max()/dx)*dx
ymin = np.ceil(y.min()/dy)*dy
ymax = np.floor(y.max()/dy)*dy


# Define X and Z coordinates of interpolated mesh
xi = np.arange(xmin,xmax+dx/2,dx)
yi = np.arange(ymin,ymax+dy/2,dy)
nx = len(xi)
ny = len(yi)

# Perform linear interpolation of the data given on original positions (x,z)
# on a grid defined by (xi,zi) and plot result (velocity vs depth)
Xi, Yi = np.meshgrid(xi, yi)

v = griddata((x,y),vel, (Xi, Yi), method='linear')
v = np.flip(v,axis=0)
fig1, ax1 = plt.subplots(1,1,figsize=(10,15))
im1 = ax1.imshow(v, cmap="rainbow", aspect='equal',\
                extent=[xmin-dx/2, xmax+dx/2, ymin-dy/2,ymax+dy/2])
ax1.set_title(f"velocities from {dir0}")
ax1.set_xlabel("distance [m]")
ax1.set_ylabel("depth [m]")
divider = make_axes_locatable(ax1)
cax1 = divider.append_axes("right", size="2%", pad=0.05)
cbar1 = plt.colorbar(im1,shrink=0.9, cax=cax1)
cbar1.set_label("velocity [m/s]")

# Calculate TWTs on grid and plot results
swt = np.zeros_like(v)
swt[1:,:] = np.nan
vv = swt.copy()
for i in range(nx):
    jt = -1
    jv = 0
    for j in range(ny-1):
        if np.isnan(v[j,i]):
            jv += 1
        else:
            break
    for j in range(jv,ny-1):
        if np.isnan(v[j,i]): break
        jt += 1
        swt[jt+1,i] = swt[jt,i]+dy/v[j,i]
        vv[jt,i] = v[j,i]
twt = swt*2.
fig2, ax2 = plt.subplots(1,1,figsize=(10,15))
im2 = ax2.imshow(twt*1000, cmap="rainbow", aspect='equal',\
                extent=[xmin-dx/2, xmax+dx/2, ymin-dy/2,ymax+dy/2])
ax2.set_title(f"TWT from {dir0}")
ax2.set_xlabel("distance [m]")
ax2.set_ylabel("depth [m]")
divider = make_axes_locatable(ax2)
cax2 = divider.append_axes("right", size="2%", pad=0.05)
cbar2 = plt.colorbar(im2,shrink=0.9, cax=cax2)
cbar2.set_label("TWT [ms]")

# Calculate RMS velocities on grid and plot results
v2 = vv**2
dswt = swt[1:,:]-swt[:-1,:]
rms = np.zeros_like(v)
rms[:,:] = np.nan
for i in range(nx):
    vsum = 0.
    for j in range(ny-1):
        if np.isnan(swt[j+1,i]): break
        vsum += v2[j,i]*dswt[j,i]
        rms[j,i] = np.sqrt(vsum/swt[j+1,i])
fig3, ax3 = plt.subplots(1,1,figsize=(10,15))
im3 = ax3.imshow(rms, cmap="rainbow", aspect='equal',\
                extent=[xmin-dx/2, xmax+dx/2, ymin-dy/2,ymax+dy/2])
ax3.set_title(f"RMS velocities from {dir0}")
ax3.set_xlabel("distance [m]")
ax3.set_ylabel("depth [m]")
divider = make_axes_locatable(ax3)
cax3 = divider.append_axes("right", size="2%", pad=0.05)
cbar3 = plt.colorbar(im3,shrink=0.9, cax=cax3)
cbar3.set_label("RMS velocity [m/s]")

x_cmp_start = x_cmp1 + (n_cmp_start-1)*dx_cmp
dx_n_cmp = d_n_cmp*dx_cmp
x_cmp_end = x_cmp1 + (n_cmp_end-1)*dx_cmp
x_cmps = np.arange(x_cmp_start,x_cmp_end+dx_n_cmp/2,dx_n_cmp)
n_cmps = np.arange(n_cmp_start,n_cmp_end+d_n_cmp/2,d_n_cmp, dtype=int)
with open("velocities.vel","w") as fo:
    fo.write("3-D: NO 0.00 2000 0\n")
    for i,n in enumerate(n_cmps):
        fo.write(f"CDP: {n} X: {x_cmps[i]:0.3f} Y: 0.000 8\n")
        for j in range(0,ny-1,n_dy):
            if np.isnan(rms[j,n]):
                break
            fo.write(f"{twt[j+1,n]*1000:0.6f} {rms[j,n]:0.6f} {vv[j,n]:0.6f} 0\n")
        j = np.where(np.isnan(rms[:,n]))[0][0]-1
        vlast = np.sqrt((rms[j,n]**2*twt[j+1,n]+vv[j,n]**2*(tmax-twt[j+1,n]))/tmax)
        fo.write(f"{tmax*1000:0.6f} {vlast:0.6f} {vv[j,n]:0.6f} 0\n")

ntimes = int(tmax/0.001)
depths = np.zeros((ntimes,nx))
yy = -np.flip(yi)
dt = 0.001
for j in range(nx):
    times = twt[:,j]
    nt = len(times)-1
    for i in range(len(times)):
        if np.isnan(times[i]):
            times[i] = 0.
            continue
        break
    for k in range(nt,i,-1):
        if np.isnan(times[k]): continue
        break
    times = times[:k]
    nt = len(times)-1
    for i,t in enumerate(np.arange(dt,tmax,dt)):
        pos = np.argmin(abs(times-t))
        if times[pos]<t and pos<nt:
            pos += 1
        depths[i+1,j] = yy[pos-1]+(yy[pos]-yy[pos-1])/(times[pos]-times[pos-1])*(t-times[pos])
fig4, ax4 = plt.subplots(1,1,figsize=(15,10))
im4 = ax4.imshow(depths, cmap="rainbow", aspect='auto',\
                extent=[xmin-dx/2, xmax+dx/2, tmax+dt/2, -dt/2])
ax4.set_title(f"Depth for given TWT from {dir0}")
ax4.set_xlabel("distance [m]")
ax4.set_ylabel("TWT [s]")
divider = make_axes_locatable(ax4)
cax4 = divider.append_axes("right", size="2%", pad=0.1)
cbar4 = plt.colorbar(im4,shrink=0.9, cax=cax4)
cbar4.set_label("Depth [m]")

