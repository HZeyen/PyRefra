# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 09:31:32 2022

@author: Hermann Zeyen
         University Paris-Saclay

Program reads first a file that contains coordinates of seismic positions (shots
and receivers). This file has the following structure:
    column 1: x_coor: longitude (UTM or decimal geographical)
    column 2: y_coor: latitude (same coordinate system as longitude)
    column 3: z_coor: altitude [m]
    column 4: x_line: x-coordinate along the profile used for the corresponding
              position (usually, the first pont will have position 0 m)
Then, the program reads the file with velocity distribution from refraPy
tomography inversion (usually, this file is called "velocities.dat"). This file
has the following structure:
    column 1: x: X_position along profile (corresponding to the coordiante system
              of x_line - see first file)
    column 2: z: Depth of velocity point (negative values)
    column 3: v: velocity
The program searches the position of x in the array x_line and interpolates its
position inside x_coor and y_coor. Then it adds z to the interpolated z_coor
position and writes new coordinates and v to file "model_projected.dat"

"""

import os
import numpy as np
#import matplotlib.pyplot as plt
import tkinter as tk
from tkinter.filedialog import askopenfilename
from pyproj import Proj, transform

Pro_g = Proj(init='epsg:4326')
Pro_u = Proj(init='epsg:32631')
Pro_l = Proj(init='epsg:2154')


dir0 = r"H:/Seg2Dat/Fontaines-Salees/2022"
os.chdir(dir0)
root = tk.Tk()
#root.withdraw()
coor_file = askopenfilename(title="Coordinates file",initialdir=dir0,
                    filetypes=(("txt Files","*.txt"),("All Files","*.*")))
with open(coor_file, "r") as f:
    x_coor = []
    y_coor = []
    z_coor = []
    x_line = []
    while True:
        try:
            nums = f.readline().split()
            x_coor.append(float(nums[0]))
            y_coor.append(float(nums[1]))
            z_coor.append(float(nums[2]))
            x_line.append(float(nums[3]))
            if x_coor[-1] < 1000:
                xx, yy = transform(Pro_g,Pro_u,x_coor[-1],y_coor[-1])
                x_coor[-1] = xx
                y_coor[-1] = yy
        except:
            break
x_coor = np.array(x_coor)
y_coor = np.array(y_coor)
z_coor = np.array(z_coor)
x_line = np.array(x_line)

model_file = askopenfilename(title="Model file",initialdir=dir0,
                    filetypes=(("txt Files","*.dat"),("All Files","*.*")))
dir0 = os.path.dirname(model_file)
os.chdir(dir0)
root.destroy()

with open(model_file, "r") as f_in:
    with open("model_projected.dat", "w") as f_out:
        while True:
            try:
                nums = f_in.readline().split()
                x = float(nums[0])
                z = float(nums[1])
                v = float(nums[2])
                pos = min(np.argmin(abs(x_line-x)),len(x_line)-2)
                if x < x_line[pos]:
                    pos = max(pos-1,0)
                fac = (x-x_line[pos])/(x_line[pos+1]-x_line[pos])
                xx = x_coor[pos]+fac*(x_coor[pos+1]-x_coor[pos])
                yy = y_coor[pos]+fac*(y_coor[pos+1]-y_coor[pos])
                zz = z_coor[pos]+fac*(z_coor[pos+1]-z_coor[pos])+z
                f_out.write(f"{yy:0.3f} {xx:0.3f} {zz:0.3f} {v} {x:0.7f}\n")
            except:
                break




