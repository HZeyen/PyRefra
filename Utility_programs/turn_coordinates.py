# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 14:27:31 2024

@author: Hermann Zeyen
         University Paris-Saclay
    
    Read PyRefra geometry files shots.geo and receivers.geo and rotate the
    coordinates of both together such that the line in average runs in X-direction
    
    User should set the names of the files to be treated in file_r (line 23, receiver file)
    and file_s (line 24, shot point file).
    Also the folder should be set in dir0 (line 22).
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.linear_model import LinearRegression

dir0 = r"E:/Seg2Dat/Feroes/Eidi_21_07_23"
file_r = "receivers_GPS.geo"
file_s = "shots_GPS.geo"
os.chdir(dir0)

data_r = np.loadtxt(file_r)
data_s = np.loadtxt(file_s)

nrec = data_r.shape[0]
nsht = data_s.shape[0]

x = np.concatenate((data_r[:,1],data_s[:,1]))
y = np.concatenate((data_r[:,2],data_s[:,2]))
z = np.concatenate((data_r[:,3],data_s[:,3]))

dx = x.max()-x.min()
dy = y.max()-y.min()
z -= z.min()

if dx > dy:
    n = np.argmin(x)
    x -= x.min()
    y -= y[n]
else:
    n = np.argmin(y)
    y -= y.min()
    x -= x[n]
neg_x = False
if x[-1] < x[0]:
    neg_x = True

model = LinearRegression().fit(x.reshape(-1, 1), y.reshape(-1, 1))
slope = model.coef_[0][0]
intercept = model.intercept_[0]

ang = np.arctan(slope)
if neg_x:
    ang += np.pi
c = np.cos(ang)
s = np.sin(ang)
xx = x*c + y*s
yy = x*s - y*c

with open("receivers.geo","w") as fo:
    for i in range(nrec):
        fo.write(f"{int(data_r[i,0])} {xx[i]:0.3f} {yy[i]:0.3f} {z[i]:0.3f}\n")

with open("shots.geo","w") as fo:
    for i in range(nsht):
        j = i+nrec
        fo.write(f"{int(data_s[i,0])} {xx[j]:0.3f} {yy[j]:0.3f} {z[j]:0.3f}\n")

fig = plt.figure(figsize=(15,13),tight_layout=True)
gs = GridSpec(15, 13, figure=fig)

if dx>dy:
    ax_ori = fig.add_subplot(gs[:6, :])
    ax_turned = fig.add_subplot(gs[7:9, :])
    ax_topo = fig.add_subplot(gs[10:15, :])
else:
    ax_ori = fig.add_subplot(gs[:, :4])
    ax_turned = fig.add_subplot(gs[:7, 5:])
    ax_topo = fig.add_subplot(gs[8:, 5:])
ax_ori.plot(x[:nrec],y[:nrec],"r+")
ax_ori.plot(x[nrec:],y[nrec:],"b*")
ax_ori.plot(x[0],y[0],"g*",ms=20)
ax_ori.set_xlabel("Easting [m]")
ax_ori.set_ylabel("Northing [m]")
ax_ori.set_title("Original coordinates")
ax_ori.set_aspect('equal')


ax_turned.plot(xx[:nrec],yy[:nrec],"r+")
ax_turned.plot(xx[nrec:],yy[nrec:],"b*")
ax_turned.plot(xx[0],yy[0],"g*",ms=20)
ax_turned.set_xlabel("Easting [m]")
ax_turned.set_ylabel("Northing [m]")
ax_turned.set_title("Turned coordinates")
ax_turned.set_aspect('equal')

ax_topo.plot(xx[:nrec],z[:nrec],"r+",label="receivers")
ax_topo.plot(xx[nrec:],z[nrec:],"b*",label="shots")
ax_topo.plot(xx[0],z[0],"g*",ms=20,label="first point")
ax_topo.set_xlabel("Easting [m]")
ax_topo.set_ylabel("Topo [m]")
ax_topo.set_title("Topography")
ax_topo.legend(bbox_to_anchor=(1,1), loc="upper right")
   