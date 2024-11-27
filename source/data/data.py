import os
import sys
from copy import deepcopy
import struct
import numpy as np
from PyQt5 import QtWidgets
from obspy.io.seg2 import seg2
from obspy.io.segy.core import _read_segy
from obspy.core.stream import Stream
from obspy.io.segy.segy import SEGYTraceHeader
from obspy.io.segy.segy import SEGYBinaryFileHeader


class Data():
    """
    Contains methods linked to data input and output:
        __init__
        readData
        getFileCorrections
        getReceiverCorrections
        mute_before_pick
        tr_head_seg2_y
        saveSEGY
        saveSU
        saveSEG2
            print_float
        seg2_write
        saveBinary
        saveASCII
        saveHeader

    """

    def __init__(self, main):
        self.main = main
        self.st = []
        self.st_ori = []
        self.save_su = False
        self.time_0 = []
        self.t0 = 0.
        self.dt = 0.
        self.tmax = 0.
        self.time = np.array([])
        self.nsamp = 0
        self.file_corr_dict = {}
        self.receiver_corr_dict = {}
        self.receiver_corr_flag = False
        self.general_sign = 1.

    def readData(self, files):
        """
        Read all chosen data files

        Parameters
        ----------
        files : string
            List of file names to be read

        """
        ntr0 = 0
# Read data file
        for nf, ff in enumerate(files.names):
            try:
                if files.file_type == "seg2":
                    try:
                        self.st.append(seg2._read_seg2(ff))
                        if 'RECEIVER_STATION_NUMBER' not in\
                                self.st[-1][0].stats.seg2:
                            for itr, tr in enumerate(self.st[-1]):
                                try:
# If trace headers do not contain the keyword "RECEIVER_STATION_NUMBER", define
#    this value as "CHANNEL NUMBER". In this case, normally also the key word
#    "SOURCE_STATION_NUMBER" is missing. It is defined as the file number
                                    tr.stats.seg2['RECEIVER_STATION_NUMBER'] =\
                                        tr.stats.seg2['CHANNEL_NUMBER']
                                    tr.stats.seg2['SOURCE_STATION_NUMBER'] =\
                                        files.numbers[nf]
# If even the keyword "CHANNEL NUMBER" is missing, RECEIVER_STATION_NUMBER" is
#    set as consecutive number of trace in the file
                                except:
                                    tr.stats.seg2['RECEIVER_STATION_NUMBER'] =\
                                        itr
                                    tr.stats.seg2['SOURCE_STATION_NUMBER'] =\
                                        files.numbers[nf]
                    except:
                        _ = QtWidgets.QMessageBox.critical(
                            None, "Error",
                            f"Error reading data file {ff}\n\nHas Obspy "
                            + " bug'NOTE' been corrected?\n"
                            + "   (see installation manual)\n\n Program stops",
                            QtWidgets.QMessageBox.Ok)
                        raise ValueError("Data file error")
                else:
                    self.st.append(_read_segy(ff, unpack_trace_headers=True))
# If SEGY files have been read, create seg2 header dictionary and integrate
#    important SEGY header information file by file and trace by trace
#    Important information are those used later like "DELAY", "UNIT_UNIQUE_ID",
#       RECEIVER_STATION_NUMBER and "SOURCE_STATION_NUMBER
                    for itr, tr in enumerate(self.st[-1]):
                        self.st[-1][itr].stats.seg2 = {}
                        self.st[-1][itr].stats.seg2['DELAY'] =\
                            self.st[-1][itr].stats.segy.trace_header.\
                            delay_recording_time/1000.
                        self.st[-1][itr].stats.seg2["UNIT_UNIQUE_ID"] = None
                        self.st[-1][itr].stats.seg2['RECEIVER_STATION_NUMBER']\
                            = self.st[-1][itr].stats.segy.\
                            trace_header.\
                            trace_number_within_the_original_field_record
                        self.st[-1][itr].stats.seg2['SOURCE_STATION_NUMBER'] =\
                            self.st[-1][itr].stats.segy.\
                            trace_header.energy_source_point_number
                self.time_0.append(self.st[-1][0].stats.seg2['DELAY'])
# get starting and end recording time
# The original program was written for SUMMIT 2 instruments. With introduction
# of SUMMIT X1, the extension has changes from sg2 to seg2 and the "DELAY" has
# changed sign. A key word unique to Summit X1 is "UNIT_UNIQUE_ID". Check the
# existance of this key word to set the correct time of the first sample with
# respect to the trigger time. It is supposed that all recordings have the
# same delay, therefore t0 is taken as the delay of the first trace of the
# first file
                if nf == 0:
                    self.t0 = float(self.st[0][0].stats.seg2['DELAY'])
                    if self.st[0][0].stats.seg2.get("UNIT_UNIQUE_ID") is not\
                            None:
                        self.t0 *= -1
                        self.general_sign = -1
# It is also supposed that all data have the same sampling rate
                    self.nsamp = self.st[0][0].stats.npts
                    self.dt = float(self.st[0][0].stats.delta)
                    self.tmax = self.dt * self.nsamp + self.t0
                    self.time = self.t0 + np.arange(self.nsamp)*self.dt
                    try:
                        self.getFileCorrections()
                        self.getReceiverCorrections()
                    except ValueError:
                        sys.exit()
# If corrections must be applied, do this now
                ifile = files.numbers[nf]
                if len(self.file_corr_dict) > 0 and ifile in\
                        self.file_corr_dict:
                    interp = self.file_corr_dict[ifile][4]
                    if interp > 1:
                        dt = float(self.st[-1][0].stats.delta)
                        fsamp = interp/dt
                        self.st[-1] = self.st[-1].resample(fsamp, window=None)
                        print("    smpling interval changed to "
                              + f"{float(self.st[-1][0].stats.delta)}")
                self.nsamp = max(self.st[-1][0].stats.npts, self.nsamp)
                self.dt = min(float(self.st[-1][0].stats.delta), self.dt)
                if len(self.file_corr_dict) > 0 and ifile in\
                        self.file_corr_dict:
                    nsht = self.file_corr_dict[ifile][0]
                    rec1 = self.file_corr_dict[ifile][1]
                    rec_step = self.file_corr_dict[ifile][2]
                    t_add = self.file_corr_dict[ifile][3]
# If trigger time is changed, shft data by the corresponding time upward or
#    downward and taper with zeros at the beginning or at the end.
#                    self.time_0[-1] += t_add
                    if abs(t_add) > 0:
                        print(f"     File {ifile}: time correction {t_add}")
                        if t_add < 0:
                            ndt = int(-t_add/self.dt)
                            for j in range(len(self.st[-1])):
                                self.st[-1][j].data[:self.nsamp-ndt] = \
                                    self.st[-1][j].data[ndt:self.nsamp]
                                self.st[-1][j].data[self.nsamp-ndt:] = 0
                        else:
                            ndt = int(t_add/self.dt)
                            for j in range(len(self.st[-1])):
                                self.st[-1][j].data[ndt:self.nsamp] = \
                                    self.st[-1][j].data[:self.nsamp-ndt]
                                self.st[-1][j].data[:ndt] = 0
# If receiver point numbers should be changed, do this now
                    if rec1 | rec_step:
                        print(f"     File {ifile}: correct receiver numbers")
                        for j in range(len(self.st[-1])):
                            self.st[-1][j].\
                                stats.seg2['RECEIVER_STATION_NUMBER'] =\
                                rec1+rec_step*j
# If shot point numbers should be changed, do this now
                    if nsht > 0:
                        for j in range(len(self.st[-1])):
                            self.st[-1][j].stats.seg2['SOURCE_STATION_NUMBER']\
                                = nsht
                            if j == 0:
                                print(f"     File {ifile}: shot point "
                                      + f"corrected, set to nr {nsht}")
                files.file_dict[nf]['traces'] =\
                    np.arange(len(self.st[-1]))+ntr0
                ntr0 = files.file_dict[nf]['traces'][-1]+1
            except ValueError:
                sys.exit()
# If an error happened reading the data, stop program,
            except Exception as exc:
                if nf == 0:
                    _ = QtWidgets.QMessageBox.critical(
                        None, "Error",
                        f"Error reading data file {ff}\n\nHas Obspy bug "
                        + "'NOTE' been corrected?\n"
                        "   (see installation manual)\n\n Program stops",
                        QtWidgets.QMessageBox.Ok)
                else:
                    _ = QtWidgets.QMessageBox.critical(
                        None, "Error",
                        f"Error reading data file {ff}\n\nProgram stops",
                        QtWidgets.QMessageBox.Ok)
                raise ImportError("Data file error") from exc
            print(f"Data set {nf} read from file {ff}; shot point: " +
                  f"{int(self.st[-1][0].stats.seg2['SOURCE_STATION_NUMBER'])}")
        print("All files read \n")
        self.time_0 = np.array(self.time_0)
        self.time = self.t0 + self.dt*np.arange(self.nsamp)
# Remove trace averages
        for s in self.st:
            for t in s:
                t.detrend('constant')
# Back up detrended data into array st_ori
        self.st_ori = deepcopy(self.st)

    def getFileCorrections(self):
        """
        Read file with file corrections (shot number, receiver numbers, timing,
             sampling rate) if it exists.
        File name must be "file_corrections.dat".

        """
        if os.path.isfile('file_corrections.dat'):
            check_dt = True
            warn_flag = True
            nl = 0
            with open('file_corrections.dat', 'r') as f:
                for line in f:
                    nums = line.split()
                    n = len(nums)
# If first line of file "file_corrections" is wrong ingore file
# If a later line is too short, it is supposed that the end of the file is
#    reached and some empty lines have been added.
                    if n < 2:
                        if nl == 0:
                            answer = QtWidgets.QMessageBox.warning(
                                None, "Warning",
                                "File 'file_corrections.dat' exists but has "
                                + "none or one column.\n\nFile is ignored\n",
                                QtWidgets.QMessageBox.Ignore |
                                QtWidgets.QMessageBox.Close,
                                QtWidgets.QMessageBox.Close)
                            del self.file_corr_dict
                        return
                    elif len(nums) < 6 and warn_flag:
                        answer = QtWidgets.QMessageBox.warning(
                            None, "Warning",
                            "File file_corrections.dat does not contain all "
                            + "expected 6 columns.\n\nMissing columns are "
                            + "ignored\nPressing 'Close' stops program",
                            QtWidgets.QMessageBox.Ok |
                            QtWidgets.QMessageBox.Close,
                            QtWidgets.QMessageBox.Close)
                        if answer == QtWidgets.QMessageBox.Close:
                            raise ValueError("File correction error: missing "
                                             + "columns")
                        warn_flag = False
                    nl += 1
                    nf = int(nums[0])
                    nsht_corr = int(nums[1])
                    nrec_start = 0
                    nrec_step = 0
                    time_add = 0.
                    interp = 1
                    if n < 3:
                        continue
                    if n > 2:
                        nrec_start = int(nums[2])
                    if n > 3:
                        nrec_step = int(nums[3])
                    if n > 4:
                        time_add = float(nums[4])
                    if n > 5:
                        interp = int(nums[5])
                    if abs(time_add) > self.tmax and check_dt:
                        answer = QtWidgets.QMessageBox.warning(
                            None, "Warning",
                            f"Time correction seems too big ({time_add:0.4f}s)"
                            + "\nIt should be given in seconds, not "
                            + "milliseconds.\nMaybe correct data in file "
                            + "'file_corrections.dat'.\n\n"
                            + "Ignore and continue or close program.",
                            QtWidgets.QMessageBox.Ignore |
                            QtWidgets.QMessageBox.Close,
                            QtWidgets.QMessageBox.Close)
                        if answer == QtWidgets.QMessageBox.Close:
                            raise ValueError("Time correction error.\n")
                        check_dt = False
                        sys.exit()
                    self.file_corr_dict[nf] = [nsht_corr, nrec_start,
                                               nrec_step, time_add, interp]

    def getReceiverCorrections(self):
        """
        Read file with receiver corrections (sign, muting) if it exists.
        File name must be "receiver_corrections.dat".

        """

# If file receiver_corrections.dat exists read the information used to cerrect
# specific receivers positions (reverse sign, eliminate trace, resample trace)
        if os.path.isfile('receiver_corrections.dat'):
            self.receiver_corr_flag = True
            nl = 0
            with open('receiver_corrections.dat', 'r') as f:
                for line in f:
                    nums = line.split()
# If first line of file "receiver_corrections" is wrong ingore file
# If a later line is too short, it is supposed that the end of the file is
#    reached and some empty lines have been added.
                    if len(nums) < 2:
                        if nl == 0:
                            _ = QtWidgets.QMessageBox.warning(
                                None, "Warning",
                                "File 'receiver_corrections.dat' exists but "
                                + "has none or one column.\n\nFile is ignored",
                                QtWidgets.QMessageBox.Ignore |
                                QtWidgets.QMessageBox.Close,
                                QtWidgets.QMessageBox.Close)
                            self.receiver_corr_flag = False
                            del self.receiver_corr_dict
                        return
                    ir = int(nums[0])-1
                    self.receiver_corr_dict[ir] = {}
                    a = float(nums[1])
                    self.receiver_corr_dict[ir]["amp"] = a
                    if np.isclose(a, 0.):
                        print(f"Receiver {ir} muted\n")
                    elif a < 0:
                        print(f"Receiver {ir} reversed sign\n")
        else:
            self.receiver_corr_flag = False

    # def saveASCII(self,v, file_out=None):
    #     """
    #     Save data of every file in ASCII format.

    #     Input:
    #         v (numpy array [n_traces,ndata]): data to be written to file
    #         file_out (str, optional)
    #         if given, data are stored in file file_out
    #         If not, the name is recNNNNN.asc where NNNNN is the number of the
    #            actual plotting window

    #     Returns
    #     -------
    #     None.

    #     """
    #     if not file_out:
    #         file_out = 'rec00000.asc'
    #     with open(file_out,'w') as fo:
    #         fo.write(f"{v.shape[1]} {v.shape[0]} "+\
    #                  "lines=data; columns=receivers\n")
    #         fo.write(f"{self.dt} {self.t0} dt [seconds], t0 [seconds]\n")
    #         np.savetxt(fo, np.transpose(v))

    def mute_before_pick(self, tr, t, delay):
        """
        Function mutes data of a trace measured before a time t_cut calculated
        as follows:
        If a time-break has been measured, the function searches for the
        minimum trace value in the time range pick-time+delay to
        pick-time+2*delay.
        From this minimum on backwards, a linear slope of length delay/4 is
        applied to the data to bring them gradually towards zero. All data
        before time-of-minimum minus slope-length are set to zero. The returned
        time t_cut is the time of centre of slope.
        If no time break has been picked, the trace is returned unchanged and
        t_cut is returned with zero velue.

        Parameters
        ----------
        tr : Numpy float array
            Contains data of trace
        t : float
            Picked travel time [s] with respect to the first sample.
            If t is negative, it is supposed that no time break exists
        delay : float
            Supposed to be the approximate period of the signal near its
            beginning [s].

        Returns
        -------
        tr : numpy float array
            Contains modified data of trace
        t_cut : float
            Contains cutting time [s] (see above in general comment for
            defintion)

        """
        if t < 0.:
            t_cut = 0.
        else:
            idel = int(delay/self.main.data.dt)
            i_start = int(t/self.main.data.dt)+idel
            imax = i_start+idel
            imin = np.argmin(tr[imax:imax+idel])+imax
            n_slope = int(idel/4)
            fac = np.arange(n_slope)/n_slope
            i1 = imin-n_slope
            tr[i1:imin] *= fac
            tr[:i1] = 0.
            t_cut = (i1+n_slope/2)*self.main.data.dt
        return tr, t_cut

    def tr_head_seg2_y(self, trace, tr_in_file, tr_in_shot, tr, x_fact=1,
                       t_fact=100, lag=0):
        """
        Fill SEGY trace header with information from SEG2 header
        Needs obspy.io.segy.segy.SEGYTraceHeader

        Input:
        tr_in_file: int
                    number of trace in the SEGY file to be written starting
                    with 0
        tr_in_shot: int
                    number of trace within shot gather. This value is written
                    tosegy header variable
                    trace_number_within_the_original_field_record
        trace: number of trace from list of all recorded traces starting with 0
        tr (Stream from obspy): Data of one trace
        x_fact : int, optional. Default: 1
                    factor with which to multiply distances (potences of 10)
        t_fact : int, optional. Default: 100
                    factor with which to multiply topography (potences of 10)
        lag : int
                   Time of the first sample before the time break in ms.
                   It is stored in trace header field LagA

        Output:
        tr (Stream from obspy): trace with SEGY header added
        """

        shot_nr = self.main.traces.shot[trace]
        receiver_nr = self.main.traces.receiver[trace]
        file_nr = self.main.traces.file[trace]

# Define header entries
        tr.stats.location = str(tr.stats.seg2.SOURCE_STATION_NUMBER)
        tr.stats.station = str(tr.stats.seg2.RECEIVER_STATION_NUMBER)
        offset = self.main.traces.offset[trace]
        if not hasattr(tr.stats, 'segy.trace_header'):
            tr.stats.segy = {}
            tr.stats.segy.trace_header = SEGYTraceHeader()
        tr.stats.segy.trace_header.ensemble_number =\
            int(self.main.traces.ncdp[trace])+1
        tr.stats.segy.\
            trace_header.x_coordinate_of_ensemble_position_of_this_trace = \
            round(self.main.traces.xcdp[trace]*x_fact)
        tr.stats.segy.trace_header.\
            y_coordinate_of_ensemble_position_of_this_trace = 0
        tr.stats.segy.trace_header.\
            distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group\
            = round(offset*x_fact)
# SEGY supports maximum 32767 samples per trace (I*2). If SEG2 trace is longer
#      than that cut data at sample 32767
        if tr.stats.npts > 32767:
            tr.stats.segy.trace_header.number_of_samples_in_this_trace = 32767
            dt = 32766*tr.stats.delta
            tr = tr.slice(starttime=tr.stats.starttime,
                          endtime=tr.stats.starttime+dt)
        else:
            tr.stats.segy.trace_header.number_of_samples_in_this_trace = \
                tr.stats.npts
        tr.stats.segy.trace_header.sample_interval_in_ms_for_this_trace = \
            int(tr.stats.delta*1000000.)
# SUMMIT2 and SUMMIT_X1 have different signs for delay of first sample
# The SEG2 header word UNIT_UNIQUE_ID characterizes SUMMIT_X1 recording
# However, since we use only pretrigger data option, lagA in SU must always be
# positive, therefore there is finally no difference.
        tr.stats.segy.trace_header.lag_time_A = lag
# Anti alias filter parameters are only indicated in seg2 header when SUMMIT2
# system is used
        try:
            tr.stats.segy.trace_header.alias_filter_frequency = \
                int(tr.stats.seg2.ALIAS_FILTER[0])
        except:
            pass
        tr.stats.segy.trace_header.year_data_recorded = tr.stats.starttime.year
        tr.stats.segy.trace_header.day_of_year = tr.stats.starttime.julday
        tr.stats.segy.trace_header.hour_of_day = tr.stats.starttime.hour
        tr.stats.segy.trace_header.minute_of_hour = tr.stats.starttime.minute
        tr.stats.segy.trace_header.second_of_minute = \
            tr.stats.starttime.second
        tr.stats.segy.trace_header.energy_source_point_number = \
            self.main.traces.unique_s_pos[tr_in_file]+1
        tr.stats.segy.trace_header.\
            number_of_vertically_summed_traces_yielding_this_trace = \
            int(tr.stats.seg2.STACK)
        tr.stats.segy.trace_header.original_field_record_number = file_nr+1
        tr.stats.segy.trace_header.shotpoint_number = \
            int(tr.stats.seg2.SOURCE_STATION_NUMBER)
        tr.stats.segy.trace_header.\
            trace_number_within_the_original_field_record = tr_in_shot
        tr.stats.segy.\
            trace_header.geophone_group_number_of_roll_switch_position_one\
            = self.main.traces.unique_r_pos[tr_in_file]+1
        tr.stats.segy.trace_header.trace_sequence_number_within_line =\
            trace+1
        tr.stats.segy.trace_header.trace_sequence_number_within_segy_file =\
            tr_in_file+1
        tr.stats.segy.trace_header.group_coordinate_x = \
            round(self.main.geo.rec_dict[receiver_nr]["x"]*x_fact)
        tr.stats.segy.trace_header.group_coordinate_y = \
            round(self.main.geo.rec_dict[receiver_nr]["y"]*x_fact)
        tr.stats.segy.trace_header.receiver_group_elevation = \
            -round(self.main.geo.rec_dict[receiver_nr]["z"]*t_fact)
        tr.stats.segy.trace_header.source_coordinate_x = \
            round(self.main.geo.sht_dict[shot_nr]["x"]*x_fact)
        tr.stats.segy.trace_header.source_coordinate_y = \
            round(self.main.geo.sht_dict[shot_nr]["y"]*x_fact)
        tr.stats.segy.trace_header.surface_elevation_at_source = \
            -round(self.main.geo.sht_dict[shot_nr]["z"]*t_fact)
        if x_fact > 1:
            tr.stats.segy.trace_header.\
                scalar_to_be_applied_to_all_coordinates = -x_fact
        else:
            tr.stats.segy.trace_header.\
                scalar_to_be_applied_to_all_coordinates = round(1/x_fact)
        if t_fact > 1:
            tr.stats.segy.trace_header.\
                scalar_to_be_applied_to_all_elevations_and_depths = -t_fact
        else:
            tr.stats.segy.trace_header.\
                scalar_to_be_applied_to_all_elevations_and_depths = \
                round(1/t_fact)
        tr.stats.segy.trace_header.instrument_gain_constant = \
            int(tr.stats.seg2.FIXED_GAIN)
        tr.stats.segy.trace_header.gain_type_of_field_instruments = 1
        tr.stats.segy.trace_header.coordinate_units = 1
        return tr

    def saveSEGY(self):
        """
        Function saves data of all shot points or only one into one or several
        files in SEGY or SU format.
        The data stored are the ones actually on the screen, including all
        filters and mutes.

        Needs:
            obspy.io.segy.segy.SEGYTraceHeader
            obspy.io.segy.segy.SEGYBinaryFileHeader

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        self.main.function = "save_SEGY"

# Open dialog window for storing parameters:
# You may store all data, all data except for the ones stored before the
# trigger or just the ones of the zoom presented on the screen.
# You may store only the actual shot or all shots
# If all shots are stored, you may put them all into one file called
# rec00000.sgy or each shot into its own file called recnnnnn.sgy, nnnnn being
# the shotpoint number
# Since coordinates are stored as integers, you may multiply them, e.g., with
# 100 if the precision should be cm.
# You may choose applying the last used frequency filter or not
# If there are several traces from the same shot and receiver points, you may
# store all multiple traces or only the first one found.
        results, okButton = self.main.dialog(
            ["Start_time",
             ["All data", "Start at time break", "save window",
              "Zero and mute before pick"],
             "All data in one file (y/n)",
             "Only this shot (y/n)",
             "Multiplicator for distances",
             "Multiplicator for topography",
             "Store multiple shot-receievers"],
            ["l", "r", "e", "e", "e", "e", "e"],
            ["b", 2, "y", "n", "1", "100", "n"], "Save SEGY/SU format")

        if okButton is False:
            print("SEGY saving cancelled")
            self.main.function = "main"
            return
        S_time = int(results[1])
        data_flag = results[2].lower()
        one_file_flag = data_flag == "y"
        file_flag = results[3].lower()
        all_shots_flag = file_flag == "n"
        mult_x = int(results[4])
        mult_topo = int(results[5])
        skip_multiple = results[6].lower() == "n"
        if skip_multiple:
            ismax = max(self.main.traces.shot)+1
            irmax = max(self.main.traces.receiver)+1
            stored = np.zeros((ismax, irmax))
#        filter_flag = False
        if S_time == 3:
            results, okButton = self.main.dialog(
                ["Approx. signal period [ms]"], ["e"], ["5"],
                "Length of source signal")
            delay = float(results[0])/1000.
        else:
            delay = 0.
        if S_time in (1, 3):
            lag = 0
        elif S_time == 0:
            lag = -int(self.t0*1000)
        else:
            lag = -int(self.main.window.time_plt_min)
        i_store = -1
        if all_shots_flag:
            if one_file_flag:
                nfiles = 1
            else:
                nfiles = len(self.main.geo.sht_dict)
            file_out = "rec00000"
            t_save = np.arange(self.main.traces.number_of_traces)
            it_save = 0
            nt_save = len(t_save)
            progressBar = QtWidgets.QProgressBar(self.main.window)
            self.main.window.mplvl.addWidget(progressBar)
            progressBar.show()
            progressBar.setValue(0)
        else:
            t_save = np.array(self.main.window.actual_traces, dtype=int)
            nfiles = 1
            data_flag = "y"
            one_file_flag = True
            file_out = f"rec{self.main.window.fig_plotted+1:0>5d}"
        shot_t_nr = np.zeros(len(self.main.traces.sht_pt_dict), dtype=int)

        if self.save_su:
            file_out = file_out+'.su'
        else:
            file_out = file_out+'.sgy'

# If data should be filtered, make first a back-up of the actual data, since
# the filter is always applied on array self.main.window.v
# May 2023: The filtering is now obsolete, since filtering may now be done for
#           all traces in the filter module. However, for the moment, I left
#           the copy of the data for simplicity (too many changes to be done)
        stw = self.st.copy()
# Start loop over all shot points
        it_save = 0
        for i in range(nfiles):
            print(f'Write file {i+1}: {file_out}')
            sst = Stream()
            if all_shots_flag:
                if one_file_flag:
                    tr_save = t_save
                else:
                    tr_save = t_save[self.main.traces.shot == i]
                    file_out = f"{file_out[:3]}{i+1:0>5}{file_out[-4:]}"
            else:
                tr_save = t_save
            if not hasattr(stw[i].stats, 'segy'):
                stw[i].stats.segy = {}
                stw[i].stats.segy.binary_file_header = SEGYBinaryFileHeader()
                stw[0].stats.segy.binary_file_header.data_sample_format_code =\
                    1

# Loop over all traces to be saved in the actual file
            for it, t in enumerate(tr_save):
                ifile = self.main.traces.file[t]
                itrace = self.main.traces.trace[t]
                ishot = self.main.traces.shot[t]
                irec = self.main.traces.receiver[t]
                start = stw[ifile][itrace].stats.starttime
                end = stw[ifile][itrace].stats.endtime
                if (it+1) % 50 == 0:
                    print(f"Trace {it+1} written: shot {ishot}, "
                          + f"receiver {irec}")
                if skip_multiple:
                    if stored[ishot, irec] > 0:
                        continue
                    else:
                        stored[ishot, irec] += 1
                shot_t_nr[ishot] += 1
                if not hasattr(stw[ifile][itrace].stats, 'segy.trace_header'):
                    stw[ifile][itrace].stats.segy = {}
                    stw[ifile][itrace].stats.segy.trace_header =\
                        SEGYTraceHeader()
                if S_time in (1, 3):
                    tr = stw[ifile].slice(starttime=start-self.t0,
                                          endtime=end)[itrace].copy()
                    tr.stats.segy.trace_header.delay_recording_time = 0
                    if S_time == 3:
                        if self.main.traces.npick[t] == 0:
                            off = self.main.traces.offset[t]
                            noff = np.where(np.abs(self.main.traces.offset-off)
                                            < 1.5)[0]
                            to = self.main.traces.trace[noff]
                            to = to[self.main.traces.npick[to] > 0]

                            tt = -1.
                        else:
                            tt = self.main.traces.pick_times[t][0]
                        tr.data, mute_t = self.mute_before_pick(tr.data, tt,
                                                                delay)
                elif S_time == 2:
                    tr = stw[ifile].slice(
                        starttime=start+self.main.window.time_plt_min-self.t0,
                        endtime=start+self.main.window.time_plt_max-self.t0
                        )[itrace].copy()
                    tr.stats.segy.trace_header.delay_recording_time = \
                        int((self.main.window.time_plt_min-self.t0)*1000)
                else:
                    tr = stw[ifile][itrace].copy()
# Set specific SegY header words
                it_save += 1
                if np.isclose(np.std(tr.data), 0.):
                    continue
                i_store += 1
#                tr = self.tr_head_seg2_y(i_store,t,tr,mult_x,mult_topo,lag)
                tr = self.tr_head_seg2_y(i_store, t, shot_t_nr[ishot], tr,
                                         mult_x, mult_topo, lag)
                if S_time == 3:
                    if np.isclose(mute_t, 0.):
                        tr.stats.segy.trace_header.data_use = 0
                    else:
                        tr.stats.segy.trace_header.mute_time_end_time = \
                            int(mute_t*1000.)
                        tr.stats.segy.trace_header.data_use = 1
                if S_time in (1, 3):
                    tr.stats.segy.trace_header.delay_recording_time = 0
                tr.data *= self.main.traces.amplitudes[t]
                tr.data = np.require(tr.data, dtype=np.float32)
                sst.append(tr)
                if all_shots_flag:
                    completed = int((it_save)/(nt_save)*100)
                    progressBar.setValue(completed)
            if self.save_su:
                sst.write(file_out, format='SU')
            else:
                sst.write(file_out, format='SEGY')
            print(f"file {file_out} written")
        if all_shots_flag:
            progressBar.setValue(0)
            self.main.window.mplvl.removeWidget(progressBar)
            progressBar.close()
            print("All data written\n")
        self.main.window.drawNew(True)
        print(f"{i_store} traces written")

        del stw
        del sst
        self.save_su = False
        self.main.function = "main"

    def saveSU(self):
        """
        Save data in SU format.
        Not recommended, since there seems to be a bug in obspy: not all header
        entries are stored.

        Function calls Save_SEGY, since the formats are the same except for
        the file header, which is not stored in SU format

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        self.save_su = True
        self.saveSEGY()

    def saveSEG2(self):
        """
        Function saves data of all shot points or only one into one or several
        files in SEG2 format.
        The data stored are the ones actually on the screen, including all
        filters and mutes.

        Needs:
            obspy.core.Stream

        Returns
        -------
        None.

        """
#        from obspy.core import Stream
        answer = self.main.test_function()
        if not answer:
            return
        self.main.function = "save_SEG2"
# Open dialog window for storing parameters:
#   You may store only the actual shot/file or all shots/files
#   If all shots are stored, you may put them all into one file called
#   shot_00000.seg2 or each shot into its own file called shot_nnnnn.sg2,
#   nnnnn being the shotpoint or file number
        folder = os.path.join(".", "seg2_save")
        results, okButton = self.main.dialog(
            ["Output folder", "Only this shot (y/n)", "Store by:",
             [" file", " shot"]], ["e", "e", "l", "r"],
            [folder, "n", "None", "1"], "Save SEG2 format")

        if okButton is False:
            print("SEG2 saving cancelled")
            self.main.function = "main"
            return
        folder = results[0]
        if not os.path.exists(folder):
            try:
                os.makedirs(folder)
            except:
                _ = QtWidgets.QMessageBox.warning(
                    None, "Warning",
                    f"Given output folder\n{folder}\n"
                    + "does not exists and cannot be created\n"
                    + "Try again clicking on save SEG2",
                    QtWidgets.QMessageBox.Close)
                return False

        single_shot_flag = results[1].lower() == "y"
        file_flag = int(results[3]) == 0
# If only the data gather actually on the sceen should be saved and a shot
#    gather is plotted, save this shot gather independent of the value given
#    for file_flag. If it is a file gather, save as file gather. For receiver
#    or distance gathers, it is not clear which shot or file gather number
#    should be saved, therefore give a warning message and leave function
        if single_shot_flag:
            if self.main.window.sg_flag:
                file_flag = False
                file = os.path.join(
                    folder, f"shot_{self.main.window.fig_plotted+1:0>5}.seg2")
                stream = Stream()
                sh = self.main.window.fig_plotted
                for i in range(len(self.main.traces.sht_pt_dict[sh]["trace"])):
                    ifile = self.main.traces.sht_pt_dict[sh]["file"][i]
                    irec = self.main.traces.sht_pt_dict[sh]["receiver"][i]
                    stream.append(self.st[ifile][irec])
                    if i == 0:
                        stream.stats = self.st[ifile].stats
            elif self.main.window.fg_flag:
                file_flag = True
                ifile = self.main.window.actual_shot-1
                file = os.path.join(folder, f"file_{ifile+1:0>5}.seg2")
                stream = self.st[ifile].copy()
            else:
                _ = QtWidgets.QMessageBox.warning(
                    None, "Warning",
                    "Only shot or file gathers may be saved in SEG2 format.\n"
                    + "Save all gathers or change plot to shot or file gather"
                    + "\n    before trying again", QtWidgets.QMessageBox.Close)
                self.main.function = "main"
                return False
            self.seg2_write(stream, file)
# If all files or shots should be saved, do this here
# Write all files to seg2 format with new headers
        else:
            if file_flag:
                for i, stream in enumerate(self.st):
                    file = os.path.join(folder, f"file_{i+1:0>5}.seg2")
                    self.seg2_write(stream, file)
# Write all shot gathers to seg2 format
            else:
                for sh in self.main.traces.sht_pt_dict:
                    stream = Stream()
                    file = os.path.join(folder, f"shot_{sh+1:0>5}.seg2")
                    for i in range(len(self.main.traces.
                                       sht_pt_dict[sh]["trace"])):
                        ifile = self.main.traces.sht_pt_dict[sh]["file"][i]
                        irec = self.main.traces.sht_pt_dict[sh]["receiver"][i]
                        stream.append(self.st[ifile][irec])
                        if i == 0:
                            stream.stats = self.st[ifile].stats
                    self.seg2_write(stream, file)
        self.main.function = "main"
        return True

    def seg2_write(self, st, file):
        """
        Writes a seismic stream object into a file using SEG2 format

        Parameters
        ----------
        st : Stream object
             The headers of the stream object must bedefined before passing it
             the seg2_write. Usually, this is done automatically if data had
             been read in using obspy.
             Accepted header dictionaries in st are "seg2" and "segy". If none
             of the two is present, seg2_write returns "False" and does not
             write the file.
        file : str
            Path to output file including file name.

        Returns
        -------
        True if valid header has been found and new file has been written or
        False else.

        """
        def print_float(a, max_decimal):
            """
            creates a character string from a float number with the minimum of
            decimals necessary up to a maximum number of decmals.
            e.g. 3.140000000001 will be printed = 3.14 (the 1 at the end is out
                                of the range of float precision)

            Parameters
            ----------
            a : float
                Floating point number to be written
            max_decimal : int
                Maximum numbe rof dicamals to be written

            Returns
            -------
            atxt : str
                text string of a with minimum number of necessary decimals

            """
            for i in range(max_decimal):
                if np.isclose(a, np.round(a, i)):
                    atxt = f"{np.round(a, i)}"
                    break
            return atxt

        nsamp = int(st[0].stats.npts)
        dt = float(st[0].stats["delta"])
        ntrace = len(st)
        ASCII_file_header = []
        file_head_chars = []
        if hasattr(st.stats, 'seg2'):
            fh = st.stats.seg2
            for key in fh:
                ASCII_file_header.append(f"{key} {fh[key]}")
                file_head_chars.append(len(ASCII_file_header[-1]))
        elif hasattr(st.stats, 'segy'):
            dt = st[0].stats['starttime']
            datum = dt.date.isoformat().replace('-', '/')
            ASCII_file_header.append(f"ACQUISITION_DATE {datum}/-")
            file_head_chars.append(len(ASCII_file_header[-1]))
            time = dt.time.isoformat()
            ASCII_file_header.append(f"ACQUISITION_TIME {time}/-")
            file_head_chars.append(len(ASCII_file_header[-1]))
            ASCII_file_header.append("CLIENT ")
            file_head_chars.append(len(ASCII_file_header[-1]))
            ASCII_file_header.append("COMPANY ")
            file_head_chars.append(len(ASCII_file_header[-1]))
            ASCII_file_header.append("INSTRUMENT ")
            file_head_chars.append(len(ASCII_file_header[-1]))
            ASCII_file_header.append("OBSERVER ")
            file_head_chars.append(len(ASCII_file_header[-1]))
            ASCII_file_header.append("TRACE_SORT COMMON_SOURCE")
            file_head_chars.append(len(ASCII_file_header[-1]))
            ASCII_file_header.append("UNITS METER")
            file_head_chars.append(len(ASCII_file_header[-1]))
            ASCII_file_header.append("NOTE ")
            file_head_chars.append(len(ASCII_file_header[-1]))
        file_head_len = np.sum(np.array(file_head_chars) + 3)
        head_len = 32 + len(st)*4 + file_head_len + 2
        ASCII_trace_headers = []
        nchars = []
        block_len = []
        trace_pointer = []
        for i in range(ntrace):
            ASCII_trace_headers.append([])
            nchars.append([])
            if hasattr(st[i].stats, "seg2"):
                th = st[i].stats["seg2"]
                for key in th:
                    key_flag = False
                    for k in fh:
                        if hasattr(k, key):
                            key_flag = True
                            break
                    if key_flag:
                        continue
                    if key == "CHANNEL_NUMBER":
                        ASCII_trace_headers[-1].append(f"{key} {i}")
                    elif key == "RECEIVER_LOCATION":
                        ASCII_trace_headers[-1].\
                            append(f"{key} {float(th[key]):0.3f}")
                    elif key == "RECEIVER_STATION_NUMBER":
                        ASCII_trace_headers[-1].append(f"{key} {int(th[key])}")
                    elif key == "SAMPLE_INTERVAL":
                        txt = print_float(th[key], 6)
                        ASCII_trace_headers[-1].append(f"{key} {txt}")
                    elif key == "SHOT_SEQUENCE_NUMBER":
                        ASCII_trace_headers[-1].append(f"{key} {int(th[key])}")
                    elif key == "SOURCE_LOCATION":
                        ASCII_trace_headers[-1].\
                            append(f"{key} {float(th[key]):0.3f}")
                    elif key == "SOURCE_STATION_NUMBER":
                        ASCII_trace_headers[-1].append(f"{key} {int(th[key])}")
                    else:
                        ASCII_trace_headers[-1].append(f"{key} {th[key]}")
                    nchars[-1].append(len(ASCII_trace_headers[-1][-1]))
            elif hasattr(st[i].stats, "segy"):
                segy_head = st[i].stats["segy"]["trace_header"]
                line = 1
                fac_topo = float(segy_head[
                    'scalar_to_be_applied_to_all_elevations_and_depths'])
                if fac_topo < 0:
                    fac_topo = -1./fac_topo
                fac_dist = float(segy_head[
                    'scalar_to_be_applied_to_all_coordinates'])
                if fac_dist < 0:
                    fac_dist = -1./fac_dist
                ASCII_trace_headers[-1].append(f"CHANNEL_NUMBER {i}")
                nchars[-1].append(len(ASCII_trace_headers[-1][-1]))
                if hasattr(segy_head, 'delay_recording_time'):
                    ASCII_trace_headers[-1].append(
                        "DELAY " +
                        f"{float(segy_head['delay_recording_time'])/1000.:0.3f}")
                else:
                    ASCII_trace_headers[-1].append("DELAY 0.0")
                nchars[-1].append(len(ASCII_trace_headers[-1][-1]))
                if hasattr(segy_head, 'instrument_gain_constant'):
                    ASCII_trace_headers[-1].append(
                        "FIXED_GAIN "
                        + f"{int(segy_head['instrument_gain_constant'])}")
                else:
                    ASCII_trace_headers[-1].append("FIXED_GAIN 1")
                nchars[-1].append(len(ASCII_trace_headers[-1][-1]))
                ASCII_trace_headers[-1].append(f"LINE_ID {line}")
                nchars[-1].append(len(ASCII_trace_headers[-1][-1]))
                ASCII_trace_headers[-1].append("POLARITY 1")
                nchars[-1].append(len(ASCII_trace_headers[-1][-1]))
                ASCII_trace_headers[-1].append(f"LINE_NUMBER {line}")
                nchars[-1].append(len(ASCII_trace_headers[-1][-1]))
                ASCII_trace_headers[-1].append(
                    "RECEIVER_LOCATION " +
                    f"{float(segy_head['group_coordinate_x']*fac_dist):0.3f}")
                nchars[-1].append(len(ASCII_trace_headers[-1][-1]))
                ASCII_trace_headers[-1].append(
                    "RECEIVER_STATION_NUMBER "
                    + f"{int(segy_head['trace_sequence_number_within_line'])}")
                nchars[-1].append(len(ASCII_trace_headers[-1][-1]))
                ASCII_trace_headers[-1].append(
                    f"SAMPLE_INTERVAL {float(st[i].stats.delta):0.6f}")
                nchars[-1].append(len(ASCII_trace_headers[-1][-1]))
                ASCII_trace_headers[-1].append(
                    "SHOT_SEQUENCE_NUMBER "
                    + "SHOT_SEQUENCE_NUMBER "
                    + f"{int(segy_head['original_field_record_number'])}")
                nchars[-1].append(len(ASCII_trace_headers[-1][-1]))
                ASCII_trace_headers[-1].append(
                    "SOURCE_LOCATION " +
                    f"{float(segy_head['source_coordinate_x'])*fac_dist:0.3f}")
                nchars[-1].append(len(ASCII_trace_headers[-1][-1]))
                ASCII_trace_headers[-1].append(
                    "SOURCE_STATION_NUMBER "
                    + f"{int(segy_head['energy_source_point_number'])}")
                nchars[-1].append(len(ASCII_trace_headers[-1][-1]))
                if hasattr(
                        segy_head,
                        'number_of_vertically_summed_traces_yielding_this_trace'):
                    ASCII_trace_headers[-1].append(
                        "STACK " +
                        f"{int(segy_head['number_of_vertically_summed_traces_yielding_this_trace'])}")
                else:
                    ASCII_trace_headers[-1].append("STACK 1")
                nchars[-1].append(len(ASCII_trace_headers[-1][-1]))
            else:
                return False
            block_len.append(np.sum(np.array(nchars[-1])+3)+32+2)
            if i == 0:
                trace_pointer.append(int(head_len))
            else:
                trace_pointer.append(int(trace_pointer[-1] +
                                         block_len[i-1]+4*nsamp))
        print(f"In routine seg2_write: Write to file {file}")
# Write file header
        with open(file, "wb") as f:
            bin_head = struct.pack('hhhhbbbbbbhhhhhhhhh', 14933, 1, 4*ntrace,
                                   ntrace, 1, 0, 32, 1, 10, 32, 0, 0, 0, 0, 0,
                                   0, 0, 0, 0)
            f.write(bin_head)
            for i in range(ntrace):
                f.write(struct.pack('I', trace_pointer[i]))
            for i, text in enumerate(ASCII_file_header):
                f.write(struct.pack("h", file_head_chars[i]+3))
                f.write(text.encode())
                f.write(struct.pack("b", 0))
            f.write(struct.pack("h", 0))
# Loop over traces
            for i in range(ntrace):
                bin_tr = struct.pack('hhIIbbhIIII', 17442, block_len[i],
                                     nsamp*4, nsamp, 4, 0, 0, 0, 0, 0, 0)
                f.write(bin_tr)
                for j, text in enumerate(ASCII_trace_headers[i]):
                    f.write(struct.pack("h", nchars[i][j]+3))
                    f.write(text.encode())
                    f.write(struct.pack("b", 0))
    # I don't know why in the next line the three last values are added, in DMT
    #   files, they are there
                f.write(struct.pack("bb", 0, 0))
                np.asarray(st[i].data, dtype=np.float32).tofile(f)
        return True

    def saveBinary(self):
        """
        Save data of every file in binary format (trace after trace) without
        header. Useful for FWI inversion program from Ecole des Mines

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        results, okButton = self.main.dialog(
            ["Start_time (a(ll)/0/w(indow))"], ["e"], [0], "Save FWI format")

        if okButton is False:
            print("FWI saving cancelled")
            return
        S_time = results[0]
        if S_time == "0":
            n1 = -int(self.t0/self.dt)
            n2 = int(self.nsamp)
        elif S_time == "w":
            n1 = int(self.main.window.time_plt_min/self.dt)
            n2 = int(self.main.window.time_plt_max/self.dt)
        else:
            n1 = 0
            n2 = self.nsamp
        print(n1, n2)
        nstart = 0
        for i, s in enumerate(self.st):
            nfil = self.main.files.numbers[i]
            nstart += 1
            file_out = f"rec{nfil:0>5}.dat"
            with open(file_out, 'w+b') as fo:
                for j in range(len(self.st[i])):
                    byte_arr = s[j].data[n1:n2]
                    nn = np.size(byte_arr)
                    binary_format = bytearray(byte_arr)
                    fo.write(binary_format)
            print(f"File {nstart}: {file_out}")
            print(f"   Number of data per trace stored in binary format: {nn}")
            print(f"   Number of traces {len(self.st[0])}")
        print("All files written in binary format\n")

    def saveASCII(self, file_out=None):
        """
        Save data of every file in ASCII format.

        Input:
            file_out (str, optional)
            if given, data are stored in file file_out
            If not, the name is recNNNNN.asc where NNNNN is the number of the
               actual plotting window

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        if not file_out:
            file_out = f"rec{self.main.window.fig_plotted+1:0>5}.asc"
        with open(file_out, 'w') as fo:
            fo.write(f"{self.main.window.v.shape[1]} "
                     + f"{self.main.window.v.shape[0]} "
                     + "lines=data; columns=receivers\n")
            fo.write(f"{self.dt} {self.t0} dt [seconds], t0 [seconds]\n")
            np.savetxt(fo, np.transpose(self.main.window.v))

    def saveHeader(self):
        """
        Writes headers of all traces to ASCII files. Every data file has its
        corresponding header file "header_nnnnn.dat". These header files are
        written to folder Headers. If this folder does not exist, it is
        created.

        Returns
        -------
        None.

        """
        answer = self.main.test_function()
        if not answer:
            return
        if not os.path.isdir("Headers"):
            os.makedirs("Headers")
        for n, s in enumerate(self.st):
            nf = self.main.files.numbers[n]
            h_file = os.path.join("Headers", f"header_{nf:0>5d}.dat")
            with open(h_file, "w") as fo:
                for itr, tr in enumerate(s):
                    fo.write(f"\ntrace {itr+1}\n")
                    fo.write(f"{tr.stats}\n")
        print("Headers written to folder Headers")
