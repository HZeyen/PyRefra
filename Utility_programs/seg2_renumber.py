# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 10:32:33 2022

@author: Hermann Zeyen
         University Paris-Saclay
"""

import os
import tkinter as tk
from tkinter.filedialog import askopenfilenames
import shutil
dir0 = r"H:/Seg2Dat/Fontaines-Salees/2022/L3_Ligne1-1"
add_num = 64
os.chdir(dir0)
dir1 = dir0+r"_Modif/"
try:
    os.mkdir(dir1)
except:
    pass
root = tk.Tk()
seg2_files = askopenfilenames(title="seg2 files",initialdir=dir0,
                    filetypes=(("seg2 Files","*.seg2"),("All Files","*.*")))
root.destroy()
for f in seg2_files:
    f2 = dir1+f.split("/")[-1]
    dot = f2.rfind(".")
    nf = int(f2[dot-5:dot])
    ff = f"{f2[:dot-5]}{nf+add_num:0>5d}{f2[dot:]}"
    print(f"{f} -> {ff}")
    shutil.copy2(f,ff)



