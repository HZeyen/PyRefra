# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 11:08:42 2024

@author: Hermann Zeyen
         University Paris-Saclay

Transform pick file in gimli *.sgt format into PyRefra *.dat format
Program reads the following files (fixed names), see PyRefra documentation for details:
    shots.geo
    receivers.geo
    picks.sgt

Output file is
    picks.dat

"""

import os
import numpy as np

#The next two lines may have to be modified:
#  SYS_PATH is the folder where all python program files are located
#  DIR0 is the folder where the data are located

DIR0 = r"E:/Seg2Dat/Campus_Route_Gif"
os.chdir(DIR0)

# Read geometry files
receivers = np.loadtxt("receivers.geo")
rec_dist = receivers[:,1]**2 + receivers[:,2]**2
shots = np.loadtxt("shots.geo")
shots_dist = shots[:,1]**2 + shots[:,2]**2

# Read sgt pick file
with open("picks.sgt","r") as fi:
    line = fi.readline()
    npos = int(line.split()[0])
    line = fi.readline()
    posx = np.zeros(npos)
    posy = np.zeros(npos)
    posz = np.zeros(npos)
    for i in range(npos):
        line = fi.readline().split()
        posx[i] = float(line[0])
        posy[i] = float(line[1])
        posz[i] = float(line[2])
    pos_dist = posx**2 + posy**2
    line = fi.readline()
    npick = int(line.split()[0])
    line = fi.readline()
    sht = np.zeros(npick, dtype=int)    
    rec = np.zeros(npick, dtype=int)    
    time = np.zeros(npick)    
    err = np.zeros(npick)    
    for i in range(npick):
        line = fi.readline().split()
        sht[i] = int(line[0])-1
        rec[i] = int(line[1])-1
        time[i] = float(line[2])
        err[i] = float(line[3])

# Calculate equivalences between shot and receiver numbering in sgt and geo files
shot_eq = np.zeros(npos, dtype = int)
rec_eq = np.zeros(npos, dtype = int)
for i,p in enumerate(pos_dist):
    n = np.argsort(abs(p-shots_dist))
    shot_eq[i] = shots[n[0],0]
    n = np.argsort(abs(p-rec_dist))
    rec_eq[i] = receivers[n[0],0]

# Store picks in PyRefra format
with open("picks.dat","w") as fo:
    for i in range(npick):
        fo.write(f"{shot_eq[sht[i]]} {rec_eq[rec[i]]} {time[i]:0.5f} "+\
                 f"{time[i]-err[i]:0.5f} {time[i]+err[i]:0.5f}\n")
