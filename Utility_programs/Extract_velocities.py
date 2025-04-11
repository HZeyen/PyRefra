
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 09:38:35 2023

@author: Hermann Zeyen
         University Paris-Saclay

Extraction of one of three possibilities:
    extraction of depth of velocity isoline
    extraction of velocities along a vertical section
    extraction of velocities along a horizontal section
Program reads file "velocities.dat" written by program PyRefra after
tomographic inversion.

User should set the folder where file velocities.dat is located

Other input is interactive
"""

import os
import sys
import numpy as np
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from PyQt5 import QtWidgets
sys_path = r"E:/Sources_2010/Python_programs"
if sys_path not in sys.path:
    sys.path.append(sys_path)
from Dialog_HZ import dialog

# Get velocity file
app = QtWidgets.QApplication(sys.argv)
files = list(QtWidgets.QFileDialog.getOpenFileName(
    None, "Select velocity file", "",
    filter="velocity (velocities.dat) ;; all (*.*)"))
file = files[0]
dir0 = os.path.dirname(file)
print(f"\nFolder: {dir0}")
# Set working folder
os.chdir(dir0)

# Read velocity file. Sign of Y coordinates is inverted since PyGimli uses
#      positive upward coordinate system
data = np.loadtxt(file)
x = data[:, 0]
y = -data[:, 1]
v = data[:, 2]

# Get user input
n2 = dir0.rfind("/")
n1 = dir0[:n2].rfind("/")
title = dir0[n1+1:n2]
labels = ["horizontal grid size [m]",
          "vertical grid size [m]",
          "Title",
          "\nWhat to be extracted:",
          ["Depth of constant velocity", "Vertical section",
           "Horizontal section"]]
types = ["e", "e", "e", "l", "r"]
values = [1., 0.5, title, None, 1]
results, okButton = dialog(labels, types, values, "Define grid and mode")
if not okButton:
    sys.exit()

dx = float(results[0])
dy = float(results[1])
title = results[2]
mode = int(results[4])

# Get coordinates or velocity to be extracted
if mode == 0:
    if v.max() > 1500.:
        v_search = 1500.
    else:
        v_search = (v.max()+v.min())*0.5
    results, okButton = dialog(["Velocity [m/s]"], ["e"], [v_search],
                               "extracted velocity")
    v_search = float(results[0])
elif mode == 1:
    x_extract = (x.max()+x.min())*0.5
    results, okButton = dialog(["Position of section [m]"], ["e"], [x_extract],
                               "extracted column")
    x_extract = float(results[0])
else:
    y_extract = (y.max()+y.min())*0.5
    results, okButton = dialog(["Position of section [m]"], ["e"], [y_extract],
                               "extracted row")
    y_extract = float(results[0])

# Extract grid limits as multiples of dx and dy. Adding/subtracting dx or dy
#    is eventually necessary to avoid extrapolation
xmn = np.round(x.min()/dx, 0)*dx
if xmn < x.min():
    xmn += dx
xmx = np.round(x.max()/dx, 0)*dx
if xmx > x.max():
    xmx -= dx
ymn = np.round(y.min()/dy, 0)*dy
if ymn < y.min():
    ymn += dy
ymx = np.round(y.max()/dy, 0)*dy
if ymx > y.max():
    ymx -= dy

# Interpolate velocities to regular grid
X_unique = np.arange(xmn, xmx+dx/2, dx)
Y_unique = np.arange(ymn, ymx+dy/2, dy)
X, Y = np.meshgrid(X_unique, Y_unique)  # 2D grid for interpolation
interp = LinearNDInterpolator(list(zip(x, y)), v)
V = interp(X, Y)

# Extract depth of velocity isoline
# Go through every column and find first occurrence of v > desired velocity
# Then calculate the precise position through linear interpolation
if mode == 0:
    z_iso = []
    x_iso = []
    for i in range(len(X_unique)):
        n = np.where(V[:, i] > v_search)[0]
        if len(n) > 0:
            j = n[0]
            x_iso.append(X_unique[i])
            z_iso.append((Y_unique[j]-Y_unique[j-1])/(V[j, i]-V[j-1, i]) *
                         (v_search-V[j-1, i])+Y_unique[j-1])
        else:
            continue
# Plot velocity model and depth of isoline
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    pl = ax.imshow(V, cmap="rainbow", extent=[xmn, xmx, ymx, ymn])
    ax.plot(x_iso, z_iso, "k")
    ax.set_xlabel("Distance[m]")
    ax.set_ylabel("Depth [m]")
    ax.set_title(f"{title} Isoline v={v_search:0.0f} m/s")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.2)
    cb = fig.colorbar(pl, cax=cax, label="velocity [m/s]")
    fig.savefig(f"isoline_{v_search:0.0f}.png")
# Save isoline coordinates to file
    with open(f"isoline_{v_search:0.0f}.txt", "w") as fo:
        fo.write(f"x, z for isoline v={v_search:0.0f} m/s from file {file}\n")
        for i in range(len(x_iso)):
            fo.write(f"{x_iso[i]:0.1f}  {z_iso[i]:0.1f}\n")
# Extract vertical section
elif mode == 1:
    ncol = int((x_extract-X_unique[0])/dx)
    v_col = V[:, ncol]
# Plot velocity model and vertical section
    fig = plt.figure(figsize=(15, 6), tight_layout=True)
    gs = GridSpec(15, 13, figure=fig)
    ax_mod = fig.add_subplot(gs[:, :9])
    ax_sec = fig.add_subplot(gs[:, 10:])
    pl = ax_mod.imshow(V, cmap="rainbow", extent=[xmn, xmx, ymx, ymn])
    ax_mod.plot([x_extract, x_extract], [ymx, ymn], "k", linewidth=3)
    ax_mod.set_xlabel("Distance[m]")
    ax_mod.set_ylabel("Depth [m]")
    ax_mod.set_title(f"{title} velocity model")
    divider = make_axes_locatable(ax_mod)
    cax = divider.append_axes("right", size="2%", pad=0.2)
    cb = fig.colorbar(pl, cax=cax, label="velocity [m/s]")
    ax_sec.plot(v_col, Y_unique, linewidth=3)
    ax_sec.set_xlabel("Velocity [m/s]")
    ax_sec.set_ylabel("Depth [m]")
    ax_sec.set_title(f"Section at x={x_extract:0.1f} m")
    ax_sec.set_ylim((ymx, ymn))
    ax_sec.grid(axis="both")
    fig.savefig(f"vertical_section_{x_extract:0.0f}.png")
# Save y-coordinate and velocity to file
    with open(f"vertical_section_{x_extract:0.0f}.txt", "w") as fo:
        fo.write(f"z, velocities at x={x_extract:0.1f} m from file {file}\n")
        for i, v in enumerate(v_col):
            if np.isfinite(v):
                fo.write(f"{Y_unique[i]:0.1f}  {v:0.0f}\n")
# Extract horizontal section
else:
    nrow = int((y_extract-Y_unique[0])/dy)
    v_row = V[nrow, :]
# Plot velocity model and horizontal section
    fig, ax = plt.subplots(2, 1, figsize=(12, 8))
    pl = ax[0].imshow(V, cmap="rainbow", extent=[xmn, xmx, ymx, ymn])
    ax[0].plot([xmn, xmx], [y_extract, y_extract], "k", linewidth=3)
    ax[0].set_xlabel("Distance[m]")
    ax[0].set_ylabel("Depth [m]")
    ax[0].set_title(f"{title} velocity model")
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes("right", size="2%", pad=0.2)
    cb = fig.colorbar(pl, cax=cax, label="velocity [m/s]")
    ax[1].plot(X_unique, v_row, linewidth=3)
    ax[1].set_xlabel("Distance [m]")
    ax[1].set_ylabel("Velocity [m/s]")
    ax[1].set_title(f"Section at z={y_extract:0.1f} m")
    ax[1].grid(axis="both")
    fig.savefig(f"horizontal_section_{y_extract:0.0f}.png")
# Save x-coordinate and velocity to file
    with open(f"horizontal_section_{y_extract:0.0f}.txt", "w") as fo:
        fo.write(f"x, velocities at z={y_extract:0.1f} m from file {file}\n")
        for i, v in enumerate(v_row):
            if np.isfinite(v_row[i]):
                fo.write(f"{X_unique[i]:0.1f}  {v_row[i]:0.0f}\n")
plt.show(fig)
