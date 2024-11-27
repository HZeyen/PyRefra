import os
from copy import deepcopy
import numpy as np
from PyQt5 import QtWidgets


class Traces():
    """
    Contains methods for trace treatment:

        __init__
        readMeasPicks
        add_pick
        storePicks
        readCalcPicks
        saveGimli
    """

    def __init__(self, main, data, geom):
        self.file = []
        self.trace = []
        self.shot = []
        self.shot_pos = []
        self.receiver = []
        self.receiver_pos = []
        self.unique_s_pos = []
        self.unique_r_pos = []
        self.component = []
        self.plotted = []
        self.nsample_trace = []
        self.dt_trace = []
        self.t0_trace = []
        self.npick = []
        self.pick_times = []
        self.pick_times_min = []
        self.pick_times_max = []
        self.offset = []
        self.off_round = []
        self.amplitudes = []
        self.nsamp_trace = []
        self.t0_trace = []
        self.dt_trace = []
        self.calc_picks = False
        self.save_su = False
        self.sht_rec_dict = {}
        self.geom = geom
        self.data = data
        self.main = main
        self.xcdp = []
# Prepare dictionaries for shot and receiver points. For every shot and
#    receiver point, three lists are prepared
        self.sht_pt_dict = {}
        for ns in geom.sht_dict.keys():
            self.sht_pt_dict[ns] = {}
            self.sht_pt_dict[ns]["file"] = []
            self.sht_pt_dict[ns]["trace"] = []
            self.sht_pt_dict[ns]["receiver"] = []
            self.sht_pt_dict[ns]["axes"] = []
        self.rec_pt_dict = {}
        for nr in geom.rec_dict.keys():
            self.rec_pt_dict[nr] = {}
            self.rec_pt_dict[nr]["file"] = []
            self.rec_pt_dict[nr]["trace"] = []
            self.rec_pt_dict[nr]["shot"] = []
            self.rec_pt_dict[nr]["axes"] = []
        self.off_pt_dict = {}
# Prepare dictionaries for shot and receiver points. For every shot and
#    receiver point, three lists are prepared
        ntr = -1
        for nf, s in enumerate(data.st):
            for nt, tr in enumerate(s):
                ntr += 1
                self.file.append(nf)
                self.trace.append(nt)
                nsht = int(tr.stats.seg2['SOURCE_STATION_NUMBER'])-1
                nrec = int(tr.stats.seg2['RECEIVER_STATION_NUMBER'])-1
                self.sht_rec_dict[(nsht, nrec)] = ntr
                self.shot.append(nsht)
                self.receiver.append(nrec)
                self.component.append(geom.rec_dict[nrec]["type"])
                self.nsample_trace.append(self.data.st[nf][nt].stats.npts)
                self.dt_trace.append(self.data.st[nf][nt].stats.delta)
                self.t0_trace.append(self.data.time_0[nf])
                self.plotted.append(False)
                self.npick.append(0)
                self.pick_times.append([])
                self.pick_times_min.append([])
                self.pick_times_max.append([])
                if nsht in self.sht_pt_dict:
                    self.sht_pt_dict[nsht]["file"].append(nf)
                    self.sht_pt_dict[nsht]["trace"].append(nt)
                    self.sht_pt_dict[nsht]["receiver"].append(nrec)
                else:
                    _ = QtWidgets.QMessageBox.critical(
                        None, "Error",
                        f"Shot {nsht} not found in sht_pt_dict\n\n"
                        + "File shots.geo may be incomplete.\n\nProgram stops",
                        QtWidgets.QMessageBox.Ok)
                    raise Exception("Shot point number missing.\n")
                if nrec in self.rec_pt_dict:
                    self.rec_pt_dict[nrec]["file"].append(nf)
                    self.rec_pt_dict[nrec]["trace"].append(nt)
                    self.rec_pt_dict[nrec]["shot"].append(nsht)
                else:
                    _ = QtWidgets.QMessageBox.critical(
                        None, "Error",
                        f"Receiver {nrec} not found in rec_pt_dict\n\n"
                        + "File receivers.geo may be incomplete.\n\n"
                        + "Program stops", QtWidgets.QMessageBox.Ok)
                    raise Exception("Receiver point number missing.\n")

# Check whether shot point from trace header exists in shot point dictionary
                if nsht in geom.sht_dict:
                    xs = geom.sht_dict[nsht]["x"]
                    ys = geom.sht_dict[nsht]["y"]
                    zs = geom.sht_dict[nsht]["z"]
                else:
                    _ = QtWidgets.QMessageBox.critical(
                        None, "Error",
                        f"Shot point number {nsht+1} not found in shots.geo\n"
                        + "\nFile shots.geo may be incomplete.\n\nProgram "
                        + "stops", QtWidgets.QMessageBox.Ok)
                    raise Exception("Shot point number missing.\n")
                self.unique_s_pos.append(self.geom.get_unique_position(
                    xs, ys, zs))
# Check whether receiver point from trace header exists in shot point
# dictionary
                if nrec in geom.rec_dict:
                    xr = geom.rec_dict[nrec]["x"]
                    yr = geom.rec_dict[nrec]["y"]
                    zr = geom.rec_dict[nrec]["z"]
                    self.shot_pos.append(xs)
                    self.receiver_pos.append(xr)
                    off = np.sqrt((xs-xr)**2+(ys-yr)**2++(zs-zr)**2)
                    self.xcdp.append(np.round(xs+xr, 0)*0.5)
                else:
                    _ = QtWidgets.QMessageBox.critical(
                        None, "Error",
                        f"Receiver point number {nrec+1} not found in " +
                        "receivers.geo\nFile receivers.geo may be incomplete."
                        + "\n\nProgram stops", QtWidgets.QMessageBox.Ok)
                    raise Exception("Receiver point number missing.\n")
                self.unique_r_pos.append(self.geom.get_unique_position(
                    xr, yr, zr))
# Calculate signed offset
                if xr < xs:
                    off = -off
                self.offset.append(off)
                self.off_round.append(round(off, 0))
                o = int(np.round(abs(off), 0))
                if o not in self.off_pt_dict:
                    self.off_pt_dict[o] = {}
                    self.off_pt_dict[o]["file"] = []
                    self.off_pt_dict[o]["trace"] = []
                    self.off_pt_dict[o]["shot"] = []
                    self.off_pt_dict[o]["receiver"] = []
                self.off_pt_dict[o]["file"].append(nf)
                self.off_pt_dict[o]["trace"].append(nt)
                self.off_pt_dict[o]["shot"].append(nsht)
                self.off_pt_dict[o]["receiver"].append(nrec)
                self.amplitudes.append(self.data.general_sign)
                if data.receiver_corr_flag:
                    if nrec in data.receiver_corr_dict:
                        self.amplitudes[-1] *=\
                            data.receiver_corr_dict[nrec]["amp"]
# Check whether all shot points from file "shots.geo" exist.
# Shot points that have not been recorded are eliminated from the dictionaries
        spt_dict = deepcopy(self.sht_pt_dict)
        for n in spt_dict:
            if len(spt_dict[n]["file"]) == 0:
                if n in self.sht_pt_dict:
                    del self.sht_pt_dict[n]
                if n in self.geom.sht_dict:
                    del self.geom.sht_dict[n]
# Check whether all receiver points from file "receivers.geo" exist.
# Receiver points that have not recorded data are eliminated from the
# dictionaries
        rpt_dict = deepcopy(self.rec_pt_dict)
        for n in rpt_dict:
            if len(rpt_dict[n]["file"]) == 0:
                if n in self.rec_pt_dict:
                    del self.rec_pt_dict[n]
                if n in self.geom.rec_dict:
                    del self.geom.rec_dict[n]
        del spt_dict, rpt_dict
        self.shot = np.array(self.shot, dtype=int)
        self.number_of_traces = len(self.shot)
        self.receiver = np.array(self.receiver, dtype=int)
        self.receiver_pos = np.array(self.receiver_pos, dtype=float)
        self.component = np.array(self.component, dtype=str)
        self.shot_pos = np.array(self.shot_pos, dtype=float)
        self.nsample_trace = np.array(self.nsample_trace, dtype=int)
        self.dt_trace = np.array(self.dt_trace, dtype=float)
        self.t0_trace = np.array(self.t0_trace, dtype=float)
        self.offset = np.array(self.offset, dtype=float)
        self.off_round = np.array(self.off_round, dtype=float)
        self.file = np.array(self.file, dtype=int)
        self.trace = np.array(self.trace, dtype=int)
        self.amplitudes = np.array(self.amplitudes, dtype=float)
        self.off_min = self.offset.min()
        self.off_max = self.offset.max()
        self.xcdp = np.array(self.xcdp)
        self.plotted = np.array(self.plotted, dtype=bool)
        self.unique_s_pos = np.array(self.unique_s_pos, dtype=int)
        self.unique_r_pos = np.array(self.unique_r_pos, dtype=int)
        cdps = np.sort(np.unique(self.xcdp))
        nc = len(self.xcdp)
        self.ncdp = np.zeros(nc, dtype=int)
        for i in range(nc):
            self.ncdp[i] = np.where(cdps == self.xcdp[i])[0][0]
        self.readMeasPicks()
        if os.path.isfile("calc_picks.dat"):
            self.readCalcPicks()
        del cdps

    def readMeasPicks(self):
        """
        Read measured picks from file picks.dat (own format)
        File contains one line per measured pick and 5 rows:
        shot_point_number, receiver_point_number, travel_time, lower and upper
        uncertainty limits

        Returns
        -------
        ishs (numpy int array) For every pick the shot point number
        ists (numpy int array) For every pick the receiver point number
        ts (numpy float array) Measured travel times
        ts_min (numpy float array) Lower uncertainty limit
        ts_max (numpy float array) Upper uncertainty limit

        """
        self.external_picks = False
        try:
            with open('picks.dat', 'r') as f:
                self.ishs = []
                self.ists = []
                self.ts = []
                self.ts_min = []
                self.ts_max = []
                while True:
                    try:
                        numbers = f.readline().split()
                        isht = int(numbers[0])-1
                        irec = int(numbers[1])-1
# If picks of read shot points are found, store them in arrays for those shot
# points
                        if (isht, irec) in self.sht_rec_dict:
                            itr = self.sht_rec_dict[(isht, irec)]
                            self.npick[itr] += 1
                            self.pick_times[itr].append(float(numbers[2]))
                            self.pick_times_min[itr].append(float(numbers[3]))
                            self.pick_times_max[itr].append(float(numbers[4]))
# If picks are found belonging to shot points not read in, store them into
# backup arrays. Thos picks will be written in function storePicks after
# picks belonging to active shot points.
                        else:
                            self.external_picks = True
                            self.ishs.append(isht)
                            self.ists.append(irec)
                            self.ts.append(float(numbers[2]))
                            self.ts_min.append(float(numbers[3]))
                            self.ts_max.append(float(numbers[4]))
                    except:
                        break
        except:
            print("\nNo pick file found")

    def add_pick(self):
        """
        Not yet used
        """
        pass

    def storePicks(self):
        """
        Store picks in file picks.dat

        Returns
        -------
        None.

        """

# if picks have been stored already in file picks.dat, read first the data of
#  this file. Then replace the read picks by new ones for all traces that have
#  been read in. The reason is that if not all existing files were chosen, it
#  is possible that existing picks from other files were not copied to the
#  Traces array and would be lost without this procedure.
        answer = self.main.test_function()
        if not answer:
            return
        try:
            nsht = int(np.max(list(self.geom.sht_dict.keys())))+1
            nrec = int(np.max(list(self.geom.rec_dict.keys())))+1
            npk = np.zeros((nsht, nrec), dtype=int)
            tpk = -np.ones((nsht, nrec, 5))
            tpkmn = np.zeros((nsht, nrec, 5))
            tpkmx = np.zeros((nsht, nrec, 5))
            print(f"\nstore {np.sum(self.npick)} picks into file picks.dat")
            for i, n in enumerate(self.npick):
                isht = self.shot[i]
                irec = self.receiver[i]
                if n == 0:
                    npk[isht, irec] = 0
                else:
                    for j in range(self.npick[i]):
                        tpk[isht, irec, j] = self.pick_times[i][j]
                        tpkmn[isht, irec, j] = self.pick_times_min[i][j]
                        tpkmx[isht, irec, j] = self.pick_times_max[i][j]
                        npk[isht, irec] = self.npick[i]
            with open("picks.dat", "w") as fh:
                for ish in range(nsht):
                    for ist in range(nrec):
                        if npk[ish, ist] > 0:
                            for ipk in range(npk[ish, ist]):
                                fh.write(f"{ish+1} {ist+1} "
                                         + f"{tpk[ish,ist,ipk]:0.5f} "
                                         + f"{tpkmn[ish,ist,ipk]:0.5f} "
                                         + f"{tpkmx[ish,ist,ipk]:0.5f}\n")
# If they exist, save picks read into backup arrays
                if self.external_picks:
                    for i, t in enumerate(self.ts):
                        fh.write(f"{self.ishs[i]+1} {self.ists[i]+1} "
                                 + f"{t:0.5f} "
                                 + f"{self.ts_min[i]:0.5f} "
                                 + f"{self.ts_max[i]:0.5f}\n")

            print("File picks.dat written")
        except:
            choice = QtWidgets.QMessageBox.warning(
                None, "Warning",
                "Error in saving picks.\nProbably contact to disk lost or\n"
                + "    file picks.dat is used by another program\n"
                + "Retry or close program\n\n"
                + "It is strongly recommended to exit and restart\n"
                + "the program even after using 'Retry'.",
                QtWidgets.QMessageBox.Retry | QtWidgets.QMessageBox.Close)
            if choice == QtWidgets.QMessageBox.Retry:
                return
            raise Exception("Pick saving error.\n")

    def readCalcPicks(self):
        """
        If file "calc_picks.dat" has been found, calculated
        travel times are plotted onto the seismogram section.
        calc_picks.dat is the standard file for travel times produced
        by the tomography algorithm integrated in this program.

        Returns
        -------
        None.

        """

# self.calc_picks is a flag set to plot calculated travel times
#    If its value is True, set it to False in order not to plot
#    calculated travel times any more
#    If it is False, change to True and plot calculated travel times into
#    all windows
        if self.calc_picks:
            self.calc_picks = False
            return
# First check whether file "calc_picks.dat" exists, which is the default
#   name and format if travel times have been generated by the integrated
#   tomography algorithm based on Pygimli.
        if os.path.isfile("calc_picks.dat"):
            self.calc_picks = True
            with open('calc_picks.dat', 'r') as f:
                dummy = np.asmatrix(f.read())
                nlines = int(np.size(dummy)/3)
                dummy = np.reshape(dummy, (nlines, 3))
                self.calc_t = np.zeros(self.number_of_traces)
                self.calc_t[:] = np.nan
                for i in range(nlines):
                    sht = int(dummy[i, 0]-1)
                    rec = int(dummy[i, 1]-1)
                    if (sht, rec) in self.sht_rec_dict:
                        trace = self.sht_rec_dict[(sht, rec)]
                        self.calc_t[trace] = float(dummy[i, 2])*0.001
                del dummy
# If no calculated picks are found give warning message
        else:
            _ = QtWidgets.QMessageBox.warning(
                None, "Warning", "File calc_picks.dat not found\n\n",
                QtWidgets.QMessageBox.Ok)
            self.calc_picks = False

    def saveGimli(self):
        """
        Store picks in file picks.sgt (format used in PyGimli, e.g. for
        seismic tomography)

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        if self.main.utilities.filtered:
            unc = 4*self.data.dt
        else:
            unc = 2*self.data.dt
        nsht = len(self.geom.sht_dict)
        nrec = len(self.geom.rec_dict)
        zmax = max(self.geom.sensors[:, 2])
        shts = list(self.geom.sht_dict.keys())
        recs = list(self.geom.rec_dict.keys())
        sensors = np.copy(self.geom.sensors)
        ncoor = len(sensors[:, 0])

        t = np.zeros((nsht, nrec))
        e = np.zeros((nsht, nrec))
        n = np.zeros((nsht, nrec), dtype=int)
        ps = np.zeros((nsht, nrec), dtype=int)
        pr = np.zeros((nsht, nrec), dtype=int)
        for i in range(nsht):
            ii = shts[i]
            for j in range(nrec):
                jj = recs[j]
                if (ii, jj) in self.sht_rec_dict:
                    ntr = self.sht_rec_dict[(ii, jj)]
                else:
                    continue
                if self.npick[ntr] > 0:
                    t[i, j] += max(self.pick_times[ntr][0], self.data.dt)
                    e[i, j] += (self.pick_times_max[ntr][0] -
                                self.pick_times_min[ntr][0])/2
                    if e[i, j] <= 0:
                        e[i, j] = unc
                    n[i, j] += 1
                    ps[i, j] = self.unique_s_pos[ntr]
                    pr[i, j] = self.unique_r_pos[ntr]

        ishot = []
        irecord = []
        tim = []
        err = []
        for i in range(nsht):
            for j in range(nrec):
                if n[i, j] > 0:
                    t[i, j] /= n[i, j]
                    e[i, j] /= n[i, j]
                    if ps[i, j] != pr[i, j]:
                        ishot.append(ps[i, j])
                        irecord.append(pr[i, j])
                        tim.append(t[i, j])
                        err.append(e[i, j])
# PyGIMLi has problems with positive Z-coordinates (it seems, rays are
# calculated, but everything above 0 is not plotted. Therefore, zmax is
# subtracted)
        sensors[:, 2] = sensors[:, 2]-zmax
        with open("picks.sgt", "w") as fo:
            ndat = len(tim)
            fo.write(f"{ncoor} # shot/geophone points\n")
# Write coordinates to file. For 2D inversion, pyGIMLi uses the y coordinate
# for depth. Therefore, the z coordiante (c[2] below) is written into the
# second column.
            fo.write("# x y z\n")
            for c in sensors:
                fo.write(f"{c[0]:0.2f} 0. {c[2]:0.2f}\n")
            fo.write(f"{ndat} # measurements\n")
            fo.write("# s g t err\n")
            for i in range(ndat):
                fo.write(f"{ishot[i]+1} {irecord[i]+1} {tim[i]:0.6f} "
                         + f"{err[i]:0.6f}\n")
        print(f"\n{ndat} picks written to file picks.sgt")
        with open("pick_dist.dat", "w") as fo:
            fo.write("Offset, time, uncertainty,midpoint: x,y,z\n")
            for i in range(ndat):
                ns = ishot[i]
                nr = irecord[i]
                o = np.sqrt((sensors[ns, 0]-sensors[nr, 0])**2 +
                            (sensors[ns, 1]-sensors[nr, 1])**2 +
                            (sensors[ns, 2]-sensors[nr, 2])**2)
                midx = (sensors[ns, 0]+sensors[nr, 0])*0.5
                midy = (sensors[ns, 1]+sensors[nr, 1])*0.5
                midz = (sensors[ns, 2]+sensors[nr, 2])*0.5
                fo.write(f"{o:0.2f} {tim[i]:0.6f} {err[i]:0.6f} {midx:0.2f} "
                         + f"{midy:0.2f} {midz:0.2f} \n")
