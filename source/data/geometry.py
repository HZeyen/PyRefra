import numpy as np
from PyQt5 import QtWidgets


class Geometry():
    """
    Methods linked to geometry information:

        __init__
        read_geo_file
        readGeom
        unique_positions
        get_unique_position

    """

    def __init__(self):
        self.rec_dict = {}
        self.sht_dict = {}
        self.pos_dict = {}
        self.sens_dict = {}
        self.types = []
        self.positions = []
        self.sensors = []
        self.x_dir = True
        self.d_x = 0.
        self.dx_geo = 0.
        self.xmin = 0.
        self.xmax = 0.
# dz_geo is distance between geophones for VSP data. Used in utilities.pseudo
        self.dz_geo = 0.

    def read_geo_file(self, filename):
        """
        Read geometry information for shots or receivers

        Parameters
        ----------
        filename : string
            Name of file to be read. Should be "receievrs.geo" or "shots.geo"

        Raises
        ------
        Exception
            Gives error message if file cannot be opened.
            Gives error message if less than 4 columns are found in a line

        Returns
        -------
        d: Dictionary with position number and coordinates; for receivers also
            component

        """
        try:
            data = np.loadtxt(filename)
        except:
            _ = QtWidgets.QMessageBox.critical(
                None, "Error",
                f"File {filename} cannot be opened or is missing.\n\n"
                + "Program stops", QtWidgets.QMessageBox.Ok)
            raise Exception(f"File {filename} cannot be opened or is "
                            + "missing.\n")

# receiver numbers in file receivers.geo start with 1, but Python neads
# numbering starting with 0. Therefore, the number read is reduced by 1 in the
# next line
        try:
            d = {}
            nl, nc = data.shape
            for i in range(nl):
                i_rec = int(data[i, 0]) - 1
                d[i_rec] = {}
                d[i_rec]["line"] = i
                d[i_rec]["x"] = data[i, 1]
                d[i_rec]["y"] = data[i, 2]
                d[i_rec]["z"] = data[i, 3]
                if nc > 4:
                    d[i_rec]["type"] = data[i, 4].upper()
                else:
                    d[i_rec]["type"] = "Z"
        except:
            if filename == "shots.geo":
                _ = QtWidgets.QMessageBox.critical(
                    None, "Error",
                    f"Error reading file {filename}, line {i+1}.\n"
                    + "Must have 4 columns: nr, X, Y, Z\n\n"
                    + "Program stops", QtWidgets.QMessageBox.Ok)
            else:
                _ = QtWidgets.QMessageBox.critical(
                    None, "Error",
                    f"Error reading file {filename}, line {i+1}.\n"
                    + "Must have 4 or 5 columns: nr, X, Y, Z [component]\n\n"
                    + "Program stops", QtWidgets.QMessageBox.Ok)
            raise Exception(f"File {filename} wrong format.\n")
        return d

    def readGeom(self):
        """
        Reads geometry files "receivers.geo" and "shots.geo"
        receivers.geo has 4 or 5 columns separated by space or tab, containing
            N, X, Y, Z [,type]
            N = number of receiver position
            X, Y, Z = coordinates of receiver position
            type (optional in each line) may be V, L, T, N, S, E or W for
                 vertical, longitudinal, transversal, N-S or E-W component
                 If type is not given, V is assumed
        shots.geo has 4 columns separated by space or tab, containing
            N, X, Y, Z
            N = number of shot point position
            X, Y, Z = coordinates of shot point position

        Output:
        x_rec, y_rec, z_rec: arrays with receiver coordinates ordered by
        number of receiver point
        t_rec : type of receiver with the following coding:
                -1: not defined
                 0: vertical component (default)
                 1: longitudinal component (parallel to line direction)
                 2: transversal component (perpendicular to line direction)
                 3: N-S component (absolute direction)
                 4: E-W component (absolute direction)
        x_sht, y_sht, z_sht: arrays with shot coordinates ordered by number of
        shot point
        d_x: distance between receiver positions
        x_dir: True if line mainly in X-direction, False if in Y-direction
        types = unique geophone types found (5th column of file receivers.geo
                if it exists, if not one value: "Z")
        """
# Read receiver geometry file receivers.geo
        self.rec_dict = self.read_geo_file("receivers.geo")

        for key in self.rec_dict:
            self.types.append(self.rec_dict[key]["type"])
        self.types = np.unique(self.types)

# Find profile direction (profiles goes mainly in X or in Y direction
        x = np.array([self.rec_dict[d]["x"] for d in self.rec_dict])
        y = np.array([self.rec_dict[d]["y"] for d in self.rec_dict])
        z = np.array([self.rec_dict[d]["z"] for d in self.rec_dict])
        self.positions = np.unique(np.vstack((x, y, z)).T, axis=0)
        for i in range(len(self.positions)):
            self.pos_dict[i] = {}
            self.pos_dict[i]["rec"] = []
            self.pos_dict[i]["comp"] = []
            for d in self.rec_dict:
                if np.isclose(self.rec_dict[d]["x"], self.positions[i, 0]) and\
                   np.isclose(self.rec_dict[d]["y"], self.positions[i, 1]) and\
                   np.isclose(self.rec_dict[d]["z"], self.positions[i, 2]):
                    self.pos_dict[i]["rec"].append(d)
                    self.pos_dict[i]["comp"].append(self.rec_dict[d]["type"])
                    self.rec_dict[d]["unique"] = i
        dx = x.max()-x.min()
        dy = y.max()-y.min()
        dz = np.abs(z[1:]-z[:-1])
        self.dz_geo = round(dz.max(), 1)
        if dx > dy:
            self.x_dir = True
            self.d_x = abs(x[1]-x[0])
            xmin_r = x.min()
            xmax_r = x.max()
            xsort = np.sort(x)
            self.dx_geo = abs(xsort[1]-xsort[0])
        else:
            self.x_dir = False
            self.d_x = abs(y[1]-y[0])
            xmin_r = y.min()
            xmax_r = y.max()
            for d in self.sht_dict:
                dum = self.sht_dict[d]["x"]
                self.sht_dict[d]["x"] = self.sht_dict[d]["y"]
                self.sht_dict[d]["y"] = dum
            for d in self.rec_dict:
                dum = self.rec_dict[d]["x"]
                self.rec_dict[d]["x"] = self.rec_dict[d]["y"]
                self.rec_dict[d]["y"] = dum
            xsort = np.sort(y)
            self.dx_geo = abs(xsort[1]-xsort[0])
        del x, y

# Read shot geometry file shots.geo
        self.sht_dict = self.read_geo_file("shots.geo")

# Calculate minimum and maximum values of receiver and shot positions
        x = np.array([self.sht_dict[d]["x"] for d in self.sht_dict])
        y = np.array([self.sht_dict[d]["y"] for d in self.sht_dict])
        if self.x_dir:
            xmin_s = x.min()
            xmax_s = x.max()
            xsort = np.sort(x)
            self.dx_geo = min(self.dx_geo, abs(xsort[1]-xsort[0]))
        else:
            xmin_s = y.min()
            xmax_s = y.max()
            for d in self.sht_dict:
                dum = self.sht_dict[d]["x"]
                self.sht_dict[d]["x"] = self.sht_dict[d]["y"]
                self.sht_dict[d]["y"] = dum
            xsort = np.sort(y)
            self.dx_geo = min(self.dx_geo, abs(xsort[1]-xsort[0]))
        self.xmin = min(xmin_r, xmin_s)
        self.xmax = max(xmax_r, xmax_s)
        self.unique_positions()

    def unique_positions(self):
        """
        Calculate and order unique positions of shots and receivers combined

        Returns
        -------
        None.

        """
        nsht = len(self.sht_dict)
        nrec = len(self.rec_dict)
        sensors_s = np.zeros((nsht, 3))
        sensors_r = np.zeros((nrec, 3))
        xs = np.array([self.sht_dict[d]["x"] for d in self.sht_dict])
        sensors_s[:, 0] = xs
        ys = np.array([self.sht_dict[d]["y"] for d in self.sht_dict])
        sensors_s[:, 1] = ys
        zs = np.array([self.sht_dict[d]["z"] for d in self.sht_dict])
        sensors_s[:, 2] = zs
        xr = np.array([self.rec_dict[d]["x"] for d in self.rec_dict])
        sensors_r[:, 0] = xr
        yr = np.array([self.rec_dict[d]["y"] for d in self.rec_dict])
        sensors_r[:, 1] = yr
        zr = np.array([self.rec_dict[d]["z"] for d in self.rec_dict])
        sensors_r[:, 2] = zr
        self.sensors = np.unique(np.concatenate((sensors_r, sensors_s)),
                                 axis=0)
        for i, s in enumerate(self.sensors):
            self.sens_dict[(s[0], s[1], s[2])] = i

    def get_unique_position(self, x, y, z):
        """
        Find number of a shot or receiver position within the unique list of
        geometry points

        Returns
        -------
        int
                number of the point in the list self.sensors (starting with 0)

        """
        return self.sens_dict[(x, y, z)]
