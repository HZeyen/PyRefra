# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 16:26:17 2024

@author: Hermann Zeyen
         University Paris-Saclay

collects picks from *.pik files and stores them in picks.dat format

"""
import sys
import os
if "E:/Sources_2010/Python_programs" not in sys.path:
    sys.path.append("E:/Sources_2010/Python_programs")
from PyQt5 import QtWidgets

app = QtWidgets.QApplication(sys.argv)
dir0 = r"E:/Seg2Dat/Rochechouart/Stage_2021/2021-09-28_Groupe2"
os.chdir(dir0)

print ("Choose pick files")
files = list(QtWidgets.QFileDialog.getOpenFileNames(None,\
            "Select PIK file", "",\
            filter="pik (*.pik) ;; all (*.*)"))
fil = files[0]
fil.sort()
with open("picks.dat","w") as fo:
    for f in fil:
        with open(f,"r") as fi:
            while True:
                line = fi.readline().split()
                if int(line[5]) == 0:
                    isht = int(line[6])
                elif int(line[5]) < 0:
                    break
                else:
                    irec = int(line[6])
                    t = float(line[3])/1000.
                    dt = float(line[4])/1000.
                    fo.write(f"{isht} {irec} {t:0.6f} {t-dt:0.6f} {t+dt:0.6f}\n")
                
