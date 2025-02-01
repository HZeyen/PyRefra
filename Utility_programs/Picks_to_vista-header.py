# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 13:37:30 2024

@author: Hermann Zeyen
         University Paris-Saclay
"""

import os

dir0 = r"E:/Seg2Dat/Fontaines-Salees/2023/Line1"
os.chdir(dir0)

with open("picks.dat","r") as fi:
    with open("out_picks.asc","w") as fo:
        fo.write("# SHOT_POINT_NO                  Col:     1 -    15 Decs:   3 Mult: 1.000000 [MEAN]\n")
        fo.write("# CHANNEL_NO                     Col:    16 -    30 Decs:   3 Mult: 1.000000 [MEAN]\n")
        fo.write("# DATA_FIRSTBREAK                Col:    31 -    45 Decs:   3 Mult: 1.000000 [MEAN]\n")
        ch = 0
        sh = 0
        for l in fi:
            ll = l.split()
            s = int(ll[0])
            c = int(ll[1])
            t = float(ll[2])*1000
            if s == sh:
                while ch < c:
                    ch += 1
                    if c == ch:
                        fo.write(f"{float(s):15.3f}{float(c):15.3f}{t:15.3f}\n")
                    else:
                        fo.write(f"{float(s):15.3f}{float(c):15.3f}{float(-1):15.3f}\n")
            else:
                sh += 1
                ch = 0
                while ch < c:
                    ch += 1
                    if c == ch:
                        fo.write(f"{float(s):15.3f}{float(c):15.3f}{t:15.3f}\n")
                    else:
                        fo.write(f"{float(s):15.3f}{float(c):15.3f}{float(-1):15.3f}\n")
                
                        
            
            