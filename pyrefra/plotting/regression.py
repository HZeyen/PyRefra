# -*- coding: utf-8 -*-
"""
Line fitting used by the picking and attenuation tools. Pure numpy/scipy, no Qt.

Contains the following functions:
    best_lines
    l1_regres
"""

import numpy as np
import scipy.stats
from scipy.optimize import minimize


def best_lines(x_points, y_points, n_lines=2, xoff1=0, Lnorm=2,
               origin=True, refra=True):
    """
    Function to calculate best fitting series of n_lines lines through
    points (x_points, y_points).
    First line is supposed to pass through the origin.
    Slopes have to decrease monotonously
    The function calculates first the unique absolute positions of
    x_points, which means that the lines are supposed to be symmetric
    around the origin. The function then searches the optimum points where
    to break slopes If more than 2 lines have to be fitted, it is possible
    to give a maximum distance for the points to be fitted for the two
    first lines.

    Parameters
    ----------
    x_points : numpy.float 1D array
        Contains X coordinates of points to be fitted
    y_points : numpy.float 1D array
        Contains Y coordinates of points to be fitted
    n_lines : int (optional, default: 2)
        Number of lines to be fitted
    xoff1 : float (optional, default: 0 (equivalent to all data))
        If given, for the fit of the first line, only points with offsets
        up to xoff1 are used.
    Lnorm : int (optional, default: 2)
        norm used for fitting (1 or 2)
    origin : bool(optional, default True)
        if True, first line must go through origin.
    refra : bool (optional, default True)
        if True, slopes are calculated for seismic refactions, i.e. slopes
            must be positive and decrease from one line to the next

    Returns
    -------
    best_off : numpy.float 1D array [n_lines+1]
        offsets of line breaks. First line goes from best_off[0] until
        best_off[1] inclusive second line goes from best_off[1] exclusive
        to best_off[2] inclusive etc.
    best_slope : numpy.float 1D array [n_lines]
        optimum slope of each line
    best_intercept : numpy.float 1D array [n_lines]
        optimum intercept of each line (first line has best_inter[0]=0)
    x_regres : numpy.float 1D array [number of unique offsets]
        unique offsets in increasing order
    y_regres : numpy.float 1D array [number of unique offsets]
        regression line y values for all unique offsets, ordered by offsets
    y_points_regres : numpy.float 1D array [number of data points]
        regression y coordinates at all positions x_points
    medians : Don't know anymore what this is...

    """

# check best regression lines crossing the obtained picks
# The direct wave and i_refra_max refractions are fitted
    if Lnorm < 1 or Lnorm > 2:
        print("Lnorm is not correctly defined, must be 1 or 2")
        print("Default value of 2 is used")
        Lnorm = 2
    n_points = len(x_points)
    best_slope = np.zeros(n_lines, dtype=float)
    best_intercept = np.zeros(n_lines, dtype=float)
    best_i = np.zeros(n_lines+1, dtype=int)
    best_off = np.zeros(n_lines+1, dtype=float)
# x_off contains the ordered offsets, i_off the corresponding trace number
# ind_c and ind_l contain the positions of the traces in the data stream st
    x_abs = np.abs(x_points)
    i_off = np.argsort(x_abs)
    x_off = x_abs[i_off]
    x_unique = np.unique(x_off)
    n_unique = len(x_unique)
# t_pick_test contains the picks ordered by offset
    y_test = y_points[i_off]
    y_regres = np.zeros(n_unique)
    x_regres = np.zeros(n_unique)
    medians = np.zeros(n_unique)

# First check for the direct wave, for this take the 15 traces with smalles
# offset and split picks into two groups (direct wave and first refracted wave)
    if n_lines > 2 and xoff1 > 0:
        max_off = xoff1
    else:
        max_off = x_unique[n_unique-1]
# xpk will contain the offsets of existing picks, ypk the corresponding times
# (t_env_der_pick may contain NaNs). i_tr is the trace number of each pick
    xpk = []
    ypk = []
    i_tr = []
# npts contains the number of picks until (and including) each position
    npts = np.zeros(n_unique+2, dtype=int)
    n_off = 0
    i = 0
    npts[0] = -1
    for n_off in range(n_unique):
        if x_unique[n_off] <= max_off:
            i_off_test = n_off
        while x_off[i] == x_unique[n_off]:
            if not np.isnan(y_test[i]):
                xpk.append(x_off[i])
                ypk.append(y_test[i])
                i_tr.append(i)
                npts[n_off] += 1
            i += 1
            if i >= n_points:
                break
        npts[n_off+1] = npts[n_off]
# create numpy arrays from lists
    xpk = np.array(xpk, dtype=float)
    ypk = np.array(ypk, dtype=float)
    yy = np.zeros_like(xpk, dtype=float)
    sig_best = 1E20
    i_best = -1

# If n_lines==1 calculate the best fitting line through all points
    if n_lines == 1 or n_unique < 4:
        if origin:
            slope1 = np.dot(xpk, ypk)/np.dot(xpk, xpk)
            intercept1 = 0.
            if Lnorm == 1:
                slope1 = l1_regres(xpk, ypk, [slope1])
        else:
            slope1, intercept1, _, _, _ = scipy.stats.linregress(xpk, ypk)
            if Lnorm == 1:
                par = l1_regres(xpk, ypk, [slope1, intercept1])
                slope1 = par[0]
                intercept1 = par[1]
        y_regres = x_off*slope1+intercept1
        x_regres = x_off
        y_points_regres = np.abs(x_points)*slope1+intercept1
        best_off[1] = np.max(x_off)
# Start loop over offsets to see for which point a line split gives the best
# fit
    else:
        for i in range(2, max(i_off_test-4, 3)):
            ne = npts[i_off_test]
# slope1 corresponds to the slope of the direct wave for the first i picks
            if origin:
                slope1 = np.dot(xpk[0:npts[i]+1], ypk[0:npts[i]+1]) /\
                         np.dot(xpk[0:npts[i]+1], xpk[0:npts[i]+1])
                if Lnorm == 1:
                    slope1 = l1_regres(xpk[0:npts[i]+1],
                                       ypk[0:npts[i]+1], [slope1])
                intercept1 = 0.
            else:
                slope1, intercept1, _, _, _ = \
                    scipy.stats.linregress(xpk[0:npts[i]+1],
                                           ypk[0:npts[i]+1])
                if Lnorm == 1:
                    par = l1_regres(
                        xpk[npts[i+1]:ne], ypk[npts[i+1]:ne],
                        [slope1, intercept1])
                    slope1 = par[0]
                    intercept1 = par[1]
# slope2 corresponds to the slope of the refacted wave for picks i+1 until 15
            slope2, intercept2, _, _, _ = \
                scipy.stats.linregress(xpk[npts[i+1]:ne],
                                       ypk[npts[i+1]:ne])
            if Lnorm == 1:
                par = l1_regres(xpk[npts[i+1]:ne], ypk[npts[i+1]:ne],
                                [slope2, intercept2])
                slope2 = par[0]
                intercept2 = par[1]
# yy contians the theoretical arrival times for the pick positions
            yy[0:npts[i]+1] = slope1*xpk[0:npts[i]+1]+intercept1
            yy[npts[i+1]:ne] = slope2*xpk[npts[i+1]:ne]+intercept2
# Calculate misfit
            diff = np.abs(ypk[:ne]-yy[:ne])
            sig = np.sqrt(np.dot(diff, diff)/len(diff))
# If misfit is smallest, save parameters
            if sig < sig_best and ((0 < slope2 < slope1
                                    and intercept2 > 0) or not refra):
                i_best = i
                s1_best = slope1
                s2_best = slope2
                int1_best = intercept1
                int2_best = intercept2
                sig_best = sig
# Store best fitting parameters and synthetic times (in regres_pick)
            if i_best < 0:
                i_best = i
                s1_best = slope1
                s2_best = slope2
                int1_best = intercept1
                int2_best = intercept2
                sig_best = sig
        best_i[1] = i_best
        best_slope[0] = s1_best
        best_intercept[0] = int1_best
        best_off[1] = x_unique[i_best]
        for i in range(best_i[1]+1):
            y_regres[i] = x_unique[i]*s1_best+int1_best
            x_regres[i] = x_unique[i]
        if n_lines > 2:
            i_best = -1
# Now calculate best fitting refractions
            for i_slope in range(1, n_lines-1):
                sig_best = 1E20
                n1 = npts[best_i[i_slope]]+1
                for n_off in range(best_i[i_slope]+3, n_unique-3):
                    n2 = npts[n_off]+1
                    slope1, intercept1, _, _, _ =\
                        scipy.stats.linregress(xpk[n1:n2], ypk[n1:n2])
                    if Lnorm == 1:
                        par = l1_regres(xpk[n1:n2], ypk[n1:n2],
                                        [slope1, intercept1])
                        slope1 = par[0]
                        intercept1 = par[1]
                    slope2, intercept2, r2, _, _ =\
                        scipy.stats.linregress(xpk[n2:], ypk[n2:])
                    if Lnorm == 1:
                        par = l1_regres(xpk[n2:], ypk[n2:],
                                        [slope1, intercept1])
                        slope2 = par[0]
                        intercept2 = par[1]
                    yy[n1:n2] = slope1*xpk[n1:n2]+intercept1
                    yy[n2:] = slope2*xpk[n2:]+intercept2
                    diff = np.abs(ypk[n1:]-yy[n1:])
                    sig = np.sqrt(np.dot(diff, diff)/len(diff))
                    if sig < sig_best and ((
                            0 < slope2 < slope1 and intercept1 > 0 and
                            intercept2 > 0) or not refra):
                        i_best = n_off
                        s1_best = slope1
                        s2_best = slope2
                        int1_best = intercept1
                        int2_best = intercept2
                        sig_best = sig
                    if i_best < 0:
                        i_best = i
                        s1_best = slope1
                        s2_best = slope2
                        int1_best = intercept1
                        int2_best = intercept2
                        sig_best = sig
                best_i[i_slope+1] = i_best
                best_off[i_slope+1] = x_unique[i_best]
                best_slope[i_slope] = s1_best
                best_slope[i_slope+1] = s2_best
                best_intercept[i_slope] = int1_best
                best_intercept[i_slope+1] = int2_best
                for i in range(best_i[i_slope], best_i[i_slope+1]+1):
                    y_regres[i] = x_unique[i]*s1_best+int1_best
            best_i[n_lines] = n_unique
            best_slope[n_lines-1] = x_unique[n_unique-1]
        elif n_lines == 2:
            best_slope[1] = s2_best
            best_intercept[1] = int2_best
        for i in range(best_i[n_lines-1], n_unique):
            y_regres[i] = x_unique[i]*s2_best+int2_best
            x_regres[i] = x_unique[i]
        y_points_regres = np.empty(n_points)
        y_points_regres.fill(np.nan)
        for i in range(n_points):
            no = np.where(x_unique == x_abs[i])
            y_points_regres[i] = y_regres[no]
        medians[0] = np.median(ypk[0:npts[1]+1])
        for i in range(1, n_unique):
            try:
                if npts[i+1] > npts[i]:
                    medians[i] = np.median(ypk[npts[i]+1:npts[i+1]+1])
            except (KeyError, IndexError, NameError):
                pass
    del x_abs
    r = np.dot((x_regres-x_regres.mean()), (y_regres-y_regres.mean())) /\
        (len(x_regres)*np.std(x_regres)*np.std(y_regres))
    r2 = r*r
    return best_off, best_slope, best_intercept, x_regres, y_regres, \
        y_points_regres, r2, medians


def l1_regres(x, y, par_ini=[0, 0]):
    """
    Function calculates linear regression line over points[x,y] using an L1
    norm

    Parameters
    ----------
    x : numpy.float 1D array
        Contains X coordinates of points to be fitted
    y : numpy.float 1D array
        Contains X coordinates of points to be fitted
    par_ini : float 1D list with one or two values optional, default: [0,0]
        Initial guess (usually obtained from a linear regression fit)
        if len(par_ini)==1, it is supposed that the given value is an
        initial guess of the slope and that the regression line passes
        through the origin if len(par_ini)==2, first parameter is slope,
        second intercept and both parameters are optimized

    Returns
    -------
    param : numpy.float 1D array , same size as par_ini
        optimized parameters
        if len(par_ini)==1, param contains only one value, the slope
        if len(par_ini)==2, param contains [slope,intercept]

    """

    def cost_function(params, X, y):
        return np.sum(np.abs(y-X.dot(params)))

    if len(par_ini) == 1:
        X = np.ones((len(x), 1))
    else:
        X = np.ones((len(x), 2))
        X[:, 0] = x
    param = minimize(cost_function, par_ini, args=(X, y))
    return param.x
