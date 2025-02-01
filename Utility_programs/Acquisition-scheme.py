# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 12:11:56 2023

@author: Hermann Zeyen

Program plots acquistion scheme. The following files are used:
    shots.geo (coordinates of shot points)
    receivers.geo (coordinates of receiver points)
        Both files have 4 columns: Point number, x, y, z coordinate
    acquisition_scheme.txt
        This file has one line per acquisition block, each with 6 numbers:
            first shot point number
            last shot point number
            first receiver point number
            last receiver point number
            first file number
            last file number
        All number are natural counting, i.e. first point has #1, not #0.
        The shot point and receiver point numbers are those stared in files
        shot.geo and receivers.geo respectively
"""

import os
import numpy as np
import matplotlib.pyplot as plt
dir0 = r"E:/Seg2Dat/Rouvres2"

os.chdir(dir0)
data = np.loadtxt("acquisition_scheme.txt")
shot1 = np.array(data[:,0], dtype=int)
shot2 = np.array(data[:,1], dtype=int)
rec1 = np.array(data[:,2], dtype=int)
rec2 = np.array(data[:,3], dtype=int)
file1 = np.array(data[:,4], dtype=int)
file2 = np.array(data[:,5], dtype=int)
nbloc = len(shot1)
shots = np.loadtxt("shots.geo")
sht_dict = {}
for i,n in enumerate(shots[:,0]):
    sht_dict[n] = shots[i,1:]
receivers = np.loadtxt("receivers.geo")
rec_dict = {}
for i,n in enumerate(receivers[:,0]):
    rec_dict[n] = receivers[i,1:]

fig1 = plt.figure(tight_layout=True)
ax1 = fig1.add_subplot()
for i in range(nbloc):
    y = i+1
    if i == 0:
        ax1.plot(sht_dict[shot1[i]][0],y,"vr",label="shots point numbers")
    else:
        ax1.plot(sht_dict[shot1[i]][0],y,"vr")
    ax1.text(sht_dict[shot1[i]][0],y+0.1,str(shot1[i]),horizontalalignment='center',\
             verticalalignment='bottom')
    x = (sht_dict[shot1[i]][0]+sht_dict[shot2[i]][0])*0.5
    ax1.text(x,y+0.1,f"files {file1[i]}-{file2[i]}",horizontalalignment='center',\
             verticalalignment='bottom')
    for j in range(shot1[i]+1,shot2[i]):
        ax1.plot(sht_dict[j][0],y,"vr")
    ax1.plot(sht_dict[shot2[i]][0],y,"vr")
    ax1.text(sht_dict[shot2[i]][0],y+0.1,str(shot2[i]),horizontalalignment='center',\
             verticalalignment='bottom')
    y -= 0.15
    if i == 0:
        ax1.plot(rec_dict[rec1[i]][0],y,"ob",label="receiver point numbers")
    else:
        ax1.plot(rec_dict[rec1[i]][0],y,"ob")
    ax1.text(rec_dict[rec1[i]][0],y-0.1,str(rec1[i]),horizontalalignment='center',\
             verticalalignment='top')
    for j in range(rec1[i]+1,rec2[i]):
        ax1.plot(rec_dict[j][0],y,"ob")
    ax1.plot(rec_dict[rec2[i]][0],y,"ob")
    ax1.text(rec_dict[rec2[i]][0],y-0.1,str(rec2[i]),horizontalalignment='center',\
             verticalalignment='top')
    ax1.legend(loc='lower right')
    ax1.set_ylim(0,nbloc+0.5)
    ax1.set_yticks(np.arange(nbloc+1))
    ax1.set_ylabel("Block number")
    ax1.set_xlabel("Position in line [m]")
fig1.savefig("acquisition_scheme.png")




