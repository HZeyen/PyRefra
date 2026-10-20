# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 09:38:35 2023

@author: Hermann Zeyen
         University Paris-Saclay

Extraction of one of three possibilities:
    extraction of depth of velocity isoline
    extraction of velocities along a vertical section
    extraction of velocities along a horizontal section
Program reads file "velocities.dat" written by program PyRefra after
tomographic inversion.

Input is interactive
"""

import os
import sys
import numpy as np
from scipy.interpolate import LinearNDInterpolator
import matplotlib
matplotlib.use("QtAgg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.colors as colors
from matplotlib.transforms import Bbox
from mpl_toolkits.axes_grid1 import make_axes_locatable
from PyQt5 import QtWidgets, QtCore
try:
    dir0 = r"E:\Seg2Dat\Normandie\line6\20260506-16.35"
    os.chdir(dir0)
except:
    pass


class Dialog(QtWidgets.QWidget):
    """
    Created on Tue Sep 13 12:50:26 2022
    @author: Hermann
    General purpose dialog box

    Input:
        parent (Class object):
            Certain values are needed from the parent class usually, the call
            will thus be:
                Dialog(self,labels,types,values,title)
        labels (list, string): labels explaining the different fields
            Explanatory text put beside the fields for data entrance.
            If values==None, label is interpreted as information/Warning text.
            For a series of Radiobuttons or Comboboxes, all corresponding
            labels are given as a list within the list. See descriptionof
            "types" for specific formats.
        types (list,string): In the order of "labels", indicates the type
            of field the label belongs to.  If labels[i] is itself  a list,
            this is considered as one single label. May have the following
            values:
            - "c": to define a check box;
            - "e": to define an editable box (lineEdit);
            - "l": To define a simple label, comment text;
            - "r": to define a series of radiobuttons (the corresponding list
              of labels is considered as one single label for the numbering of
              labels, types and values);
            - "b": to define a combobox (dropdown menu).

            For 'l', 'e' and 'b', capitalization does not matter. However for
            'c' and 'r' it does:
            For these two types, if capital letter is used, this means that the
            fields change when the corresponding button is clicked. Usually,
            one or more widgets will be added to the dialog box or removed from
            it. If this is the case, the label(s) (and type(s) and value(s)) of
            the optional fields must follow immediately the definition of the
            check box/radio button and the label must have the following form:
            @ntext. "@" is the sign for "optional", n is the value of the
            controlling button (may be several ciffers). For check box, n=0
            means that the optional field is activated if the combo box is
            deactivated, if n=1, the field is activated if the box is
            activated. For radio buttons, n is the number of the radio button
            (starting at 0) which must be checked for the optional field to be
            activated.
            Example:
                label="check", type="C", value=0
                label="@1optional1", type="e", value=100.
                label="@1optional2", type="e", value=-1.
                label="@0no entry", type="l", value=None
            In this case, if the check box is checked, two lineEdit boxes will
            be shown with labels "optional1" and "optional2" and default values
            100. and -1. If the check box is unchecked, a text line will be
            shown "no entry". It is, however, not necessary to give a field for
            both cases (checked and unchecked). Usually, the optional
            fields will only be shown after checking a box.
        values (list, string, float or int): initial values for LineEdit
            - Initial values for LineEdit fields (float, int or str)
            - None for comboboxes
            - >0 for check box if it should be checked from the beginning,
              0 else.
            - Number of radiobutton to be activated by default (natural
              numbering, starting at 1, not at 0).
            - For labels, it may be "b" (bold text), "i" (italic) or anything
              else, including None, for standard text.
        title (string): Title of the dialog window
                            Optional, default value: "Title"
    """
# Create window

    def __init__(self, labels, types, values=None, title="Title"):
        import numpy as np
        super(Dialog, self).__init__()
        nlab = len(labels)
        self.labels = labels
        self.types = types
        self.values = values
        self.ck_order = np.zeros(nlab, dtype=int)
        self.n_checked = 0
        self.Dfinish = False
        self.Dbutton = False
        self.dlabels = []
        self.dlines = []
        self.ckState = []
        self.ckb = []
        self.rbtn = []
        self.btngroup = []
        self.combo = []
        self.radio = False
        self.check = False
        nlab = len(self.types)
        for i in range(nlab):
            self.ckState.append(False)
        self.YesBtn = QtWidgets.QPushButton('Ok', self)
        self.YesBtn.move(10, 20*(nlab+3))
        self.CancelBtn = QtWidgets.QPushButton('Cancel', self)
        self.CancelBtn.move(150, 20*(nlab+3))
        self.mainLayout = QtWidgets.QGridLayout()
        self.setLayout(self.mainLayout)
        il_add = 0
        ilin = 0
        for i, t in enumerate(self.types):
            il = ilin + il_add
            if t.lower() == 'l':
                if self.values[i]:
                    if self.values[i].lower() == 'b':
                        self.dlabels.append(QtWidgets.QLabel
                                            ("<b>"+self.labels[i]+"</b>"))
                    elif self.values[i].lower() == 'i':
                        self.dlabels.append(QtWidgets.QLabel
                                            ("<i>"+self.labels[i]+"</i>"))
                    else:
                        self.dlabels.append(QtWidgets.QLabel(self.labels[i]))
                else:
                    self.dlabels.append(QtWidgets.QLabel(self.labels[i]))
                self.dlines.append(None)
                self.ckb.append(None)
                self.rbtn.append(None)
                self.btngroup.append(None)
                self.combo.append(None)
                self.mainLayout.addWidget(self.dlabels[ilin], il, 0, 1, 2)
                ilin += 1
            elif t.lower() == 'e':
                self.dlabels.append(QtWidgets.QLabel(self.labels[i]))
                self.dlines.append(QtWidgets.QLineEdit())
                self.ckb.append(None)
                self.rbtn.append(None)
                self.btngroup.append(None)
                self.combo.append(None)
                self.mainLayout.addWidget(self.dlabels[ilin], il, 0, 1, 1)
                self.mainLayout.addWidget(self.dlines[ilin], il, 1, 1, 1)
                try:
                    s = str(self.values[i])
                    self.dlines[-1].setText(s)
                except:
                    pass
                ilin += 1
            elif t.lower() == 'r':
                self.ckb.append(None)
                self.combo.append(None)
                self.rbtn.append([])
                self.btngroup.append(QtWidgets.QButtonGroup())
                rck = int(self.values[i])-1
                if rck < 0 or rck >= len(self.labels[i]):
                    rck = 0
                for ir, lab in enumerate(self.labels[i]):
                    self.dlabels.append(None)
                    self.dlines.append(None)
                    self.rbtn[i].append(QtWidgets.QRadioButton(lab))
                    self.rbtn[i][-1].setProperty("groupName", str(i))
                    self.btngroup[-1].addButton(self.rbtn[i][-1])
                    self.mainLayout.addWidget(self.rbtn[i][-1], il, 0, 1, 2)
                    if ir == rck:
                        self.rbtn[i][-1].setChecked(True)
                    else:
                        self.rbtn[i][-1].setChecked(False)
                    il += 1
                    ilin += 1
                if t == "R":
                    self.btngroup[-1].buttonClicked.connect(self.radio_checked)
            elif t.lower() == 'c':
                self.dlabels.append(None)
                self.dlines.append(None)
                self.rbtn.append(None)
                self.btngroup.append(None)
                self.combo.append(None)
                self.ckb.append(QtWidgets.QCheckBox(self))
                self.ckb[i].setText(self.labels[i])
                self.mainLayout.addWidget(self.ckb[i], il, 0, 1, 2)
                self.ckb[i].stateChanged.connect(self.checked)
                if values[i]:
                    self.ckb[i].setChecked(True)
                    self.ckState[i] = True
                    self.ck_order[i] = max(self.ck_order) + 1
                    self.check = False
                ilin += 1
            elif types[i].lower() == 'b':
                self.dlabels.append(None)
                self.dlines.append(None)
                self.ckb.append(None)
                self.rbtn.append(None)
                self.btngroup.append(None)
                self.combo.append(QtWidgets.QComboBox())
                for il, lab in enumerate(self.labels[i]):
                    self.combo[i].addItem(lab)
                ilin += 1
                self.mainLayout.addWidget(self.combo[i], ilin, 0, 1, 1)
                il_add += 1
        ilin += 2
        il = ilin + il_add
        self.mainLayout.addWidget(self.YesBtn, il, 0)
        self.mainLayout.addWidget(self.CancelBtn, il, 1)
        self.YesBtn.setDefault(True)
        self.YesBtn.clicked.connect(self.on_YesButton_clicked)
        self.CancelBtn.clicked.connect(self.on_CancelButton_clicked)
        self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
#        self.setModal(True)
        self.setWindowTitle(title)
        self.show()

    def radio_checked(self, button):
        self.radio = True
        self.checked_field = int(button.property("groupName"))
        self.on_YesButton_clicked()

    def checked(self, checked):
        """
        Actions executed if a check box has changed its state.
        if a box has been checked, the function searches the one which was
        checked using self.ckState as indicator (this variable contains the
        state of all check boxes before the click) Its click-order is stored
        and self.ckState is changed. In addition, for check boxes, the
        explanation text is changed and the order of checking is indicated.
        If a box is unchecked, this is stored in self.ckState and the order
        numbers indicated are changed if necessary.

        Parameters
        ----------
        checked : QtWidget.QCheckBox.checkState
            state of a checkbox after clicking into it

        Returns
        -------
        None.

        """
        from PyQt5.QtCore import Qt
# If a check box is checked, search the one which has been checked and
# do necessary changes
# If self.ckb has None value, the corresponding entry is not a checkbox
        if checked == Qt.Checked:
            for i, ck in enumerate(self.ckb):
                if not ck:
                    continue
# If self.ckb.checkState is checked after click, set ckState to True and do
# changes
                if ck.checkState() == Qt.Checked:
                    self.ckState[i] = True
# If checkbox nr i was not checked, increase the number of checked boxes
#    (n_checked) and store the order of checkin in self.ck_order
                    if self.ck_order[i] == 0:
                        self.n_checked += 1
                        self.ck_order[i] = self.n_checked
# if checkboxis checked, change the label, indicating the order in which it was
#    checked.
                        self.ckb[i].setText(f"{self.labels[i]} "
                                            f"({self.ck_order[i]})")
                        if self.types[i] == "C":
                            self.check = True
                            self.checked_field = i
                            self.on_YesButton_clicked()
                        # break
# If self.ckb.checkState is still unchecked, set ckState to Falsee
                else:
                    self.ckState[i] = False
# If click has unchecked a checkbox, do necessary changes
# If self.ckb has None value, the corresponding entry is not a checkbox
        else:
            for i, ck in enumerate(self.ckb):
                if not ck:
                    continue
# If self.ckb.checkState is still checked, set ckState to True
                if ck.checkState() == Qt.Checked:
                    self.ckState[i] = True
# If checkbox is no longer checked but it was (self.ckState), the unchecked box
#    is found
                else:
                    if self.ckState[i] is True:
                        self.ckState[i] = False
                        n = self.ck_order[i]
# reset ck_order to 0 (indicating also unchecked box)
                        self.ck_order[i] = 0
# Reset label to initial value (changes only for function "False_Color")
                        self.ckb[i].setText(self.labels[i])
# For all boxes that were checked later than the unchecked one, reduce their
#    checking order by 1
                        for j in range(len(self.ck_order)):
                            if self.ck_order[j] > n:
                                self.ck_order[j] -= 1
                                self.ckb[j].setText(f"{self.labels[j]} "
                                                    f"({self.ck_order[j]})")
                        self.n_checked -= 1
                        if self.types[i] == "C":
                            self.check = True
                            self.checked_field = i
                            self.on_YesButton_clicked()
                        # break
        self.show()

    def on_YesButton_clicked(self):
        from PyQt5.QtCore import Qt
        n_checked = 0
        for ck in self.ckb:
            if not ck:
                continue
            if ck.checkState() == Qt.Checked:
                n_checked += 1
        self.Dfinish = True
        self.Dbutton = True
        self.close()

    def on_CancelButton_clicked(self):
        self.Dfinish = True
        self.Dbutton = False
        self.close()


def dialog(labels, types, values, title="Title"):
    """
      Wrapper for class Dialog. Two buttons are shown for finishing: Ok and
      Cancel

      Parameters
      ----------
      labels : list of strings
          Explanatory text put beside the fields foe data entrance.
          If values==None, label is interpreted as information/Warning text
          For a series of radiobuttons or Comboboxes, all corresponding labels
          are given as a list within the list.
      types : list of strings (length as labels).
          Possible values:
              "c": to define a check box
              "e": to define an editable box (lineEdit)
              "l": To define a simple label, comment text
              'r': to define a series of radiobuttons (the corresponding list
                   of labels is considered as one single label for the
                   numbering of labels, types and values)
              "b": to define a combobox (dropdown menu)
        values (list, string, float or int): initial values for LineEdit
            - Initial values for LineEdit fields (float, int or str)
            - None for comboboxes
            - >0 for check box if it should be checked from the beginning,
              0 else.
            - Number of radiobutton to be activated by default (natural
              numbering, starting at 1, not at 0).
            - For labels, it may be "b" (bold text), "i" (italic) or anything
              else, including None, for standard text.
      title : str, default: "Title"
          Title of the dialog box.

      Returns
      -------
      results : list of bool
          Response of each data entrance field. Should be transformed to the
          needed data format (int, float...) in the calling function
          For radiobuttons, the returned value indicates the number of the
          active button (counting starts at 0). For checkboxes, the returned
          value is -1 if the box was not checked and the position at which
          the box was checked (starting at 0) if it has been checked
      Dbutton: bool
          If True, "Apply" button has been pressed to finish dialog, if False
          "Cancel" button has been pressed.

    """
    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(sys.argv)
    # Keep app alive between successive dialogs in script mode (VS Code).
    app.setQuitOnLastWindowClosed(False)
    while True:
        c_flag = -1
        r_flag = -1
        lab = []
        typ = []
        val = []
        field = []
        radio = []
        check = []
        for i, t in enumerate(types):
            icolon = 0
            if t == "C":
                lab.append(labels[i])
                typ.append(t)
                val.append(values[i])
                field.append(i)
                radio.append(None)
                check.append(values[i])
                c_flag = len(lab)-1
                r_flag = -1
            elif t == "R":
                lab.append(labels[i])
                typ.append(t)
                val.append(values[i])
                field.append(i)
                radio.append(values[i])
                check.append(None)
                r_flag = len(lab)-1
                c_flag = -1
            else:
                if r_flag > -1 and labels[i][0] == "@":
                    icolon = labels[i].index(":")
                    if icolon < 2:
                        print(f"Error: Field {i}: Colon missing")
                        _ = QtWidgets.QMessageBox.warning(
                            None, "Warning", "the radio button indicator is "
                            f"capitalized but field {i} does not start with @."
                            "\n\nCapitalization is removed.",
                            QtWidgets.QMessageBox.Close,
                            QtWidgets.QMessageBox.Close)
                        r_flag = -1
                        lab.append(labels[i])
                        typ.append(types[i])
                        val.append(values[i])
                        field.append(i)
                        radio.append(None)
                        check.append(None)
                        continue
                    index = int(labels[i][1:icolon])
                    icolon += 1
                    if val[r_flag] == index:
                        lab.append(labels[i][icolon:])
                        typ.append(types[i])
                        val.append(values[i])
                        field.append(i)
                        radio.append(None)
                        check.append(None)
                    continue
                if c_flag > -1 and labels[i][0] == "@":
                    icolon = labels[i].index(":")
                    if icolon < 2:
                        print(f"Error: Field {i}: Colon missing")
                        _ = QtWidgets.QMessageBox.warning(
                            None, "Warning", "The check box indicator is "
                            f"capitalized but field {i} does not start with @."
                            "\n\nCapitalization is removed.",
                            QtWidgets.QMessageBox.Close,
                            QtWidgets.QMessageBox.Close)
                        r_flag = -1
                        lab.append(labels[i])
                        typ.append(types[i])
                        val.append(values[i])
                        field.append(i)
                        radio.append(None)
                        check.append(None)
                        continue
                    index = int(labels[i][1:icolon])
                    if val[c_flag] > 0:
                        vv = 1
                    else:
                        vv = 0
                    icolon += 1
                    if vv == index:
                        lab.append(labels[i][icolon:])
                        typ.append(types[i])
                        val.append(values[i])
                        field.append(i)
                        radio.append(None)
                        check.append(None)
                    continue
                lab.append(labels[i][icolon:])
                typ.append(types[i])
                val.append(values[i])
                field.append(i)
                radio.append(None)
                check.append(None)
                c_flag = -1
                r_flag = -1
        D = Dialog(lab, typ, val, title)
        D.Dfinish = False
        while (D.Dfinish is not True):
            QtCore.QCoreApplication.processEvents()

        results = []
        for i in range(len(labels)):
            results.append(None)
        iline = 0
        lv = len(values)
        if lv > 0:
            for it, t in enumerate(typ):
                if t.lower() == "e":
                    results[field[it]] = D.dlines[iline].text()
                    values[field[it]] = D.dlines[iline].text()
                    iline += 1
                elif t.lower() == "r":
                    for i in range(len(lab[it])):
                        iline += 1
                        if D.rbtn[it][i].isChecked():
                            results[field[it]] = i
                            values[field[it]] = i
                elif t.lower() == 'c':
                    results[field[it]] = D.ck_order[it]-1
                    if results[field[it]] < 0:
                        values[field[it]] = 0
                    else:
                        values[field[it]] = 1
                    iline += 1
                elif t.lower() == "b":
                    results[field[it]] = D.combo[it].currentIndex()
                    iline += 1
                else:
                    iline += 1
        if D.radio is False and D.check is False:
            break
        it = D.checked_field
        if D.radio:
            D.radio = False
            values[field[it]] = results[field[it]]+1
        if D.check:
            D.check = False
            if results[field[it]] < 0:
                values[field[it]] = 0
            else:
                values[field[it]] = results[field[it]]+1
    return results, D.Dbutton


# Get velocity file
app = QtWidgets.QApplication.instance()
if app is None:
    app = QtWidgets.QApplication(sys.argv)
app.setQuitOnLastWindowClosed(False)
files = list(QtWidgets.QFileDialog.getOpenFileName(
    None, "Select velocity file", "",
    filter="velocity (velocities.dat) ;; all (*.*)"))
file = files[0]
if not file:
    sys.exit()
dir0 = os.path.dirname(file)
print(f"\nFolder: {dir0}")
# Set working folder
os.chdir(dir0)

# Read velocity file. Sign of Y coordinates is inverted since PyGimli uses
#      positive upward coordinate system
data = np.loadtxt(file)
x = data[:, 0]
y = -data[:, 1]
v = data[:, 2]
v_min = np.floor(v.min()/10.)*10.
v_max = np.ceil(v.max()/10.)*10.

# Get user input
n2 = dir0.rfind("/")
n1 = dir0[:n2].rfind("/")
title = dir0[n1+1:n2]
labels = ["horizontal grid size [m]",
          "vertical grid size [m]",
          "Title",
          "\nWhat to be extracted:",
          ["Depth of constant velocity", "Vertical section",
           "Horizontal section"],
          "Horizontal=vertical scale"]
types = ["e", "e", "e", "l", "r", "c"]
values = [1., 0.5, title, None, 1, -1]
results, okButton = dialog(labels, types, values, "Define grid and mode")
if not okButton:
    sys.exit()

dx = float(results[0])
dy = float(results[1])
title = results[2]
mode = int(results[4])
scaling = int(results[5]) > 0

# Get coordinates or velocity to be extracted
if mode == 0:
    if v.max() > 1500.:
        v_search = 1500.
    else:
        v_search = (v.max()+v.min())*0.5
    results, okButton = dialog(
        ["Velocity [m/s]",
         "Minimum velocity color scale [m/s]",
         "Maximum velocity color scale {m/s]",
         "if v_min==v_max, extreme values are used"],
        ["e", "e", "e", "l"], [v_search, v_min, v_max, None],
        "extracted velocity")
    if not okButton:
        sys.exit()
    v_search = float(results[0])
    v_min = float(results[1])
    v_max = float(results[2])
elif mode == 1:
    x_extract = np.round((x.max()+x.min())*0.5, 1)
    results, okButton = dialog(
        ["Position of section [m]",
         "Minimum velocity color scale [m/s]",
         "Maximum velocity color scale {m/s]",
         "if v_min==v_max, extreme values are used"],
        ["e", "e", "e", "l"], [x_extract, v_min, v_max, None],
        "extracted column")
    if not okButton:
        sys.exit()
    x_extract = float(results[0])
    v_min = float(results[1])
    v_max = float(results[2])
else:
    y_extract = np.round((y.max()+y.min())*0.5, 1)
    results, okButton = dialog(
        ["Position of section [m]",
         "Minimum velocity color scale [m/s]",
         "Maximum velocity color scale {m/s]",
         "if v_min==v_max, extreme values are used"],
        ["e", "e", "e", "l"], [y_extract, v_min, v_max, None],
        "extracted row")
    if not okButton:
        sys.exit()
    y_extract = float(results[0])
    v_min = float(results[1])
    v_max = float(results[2])

# Extract grid limits as multiples of dx and dy. Adding/subtracting dx or dy
#    is eventually necessary to avoid extrapolation
xmn = np.round(x.min()/dx, 0)*dx
if xmn < x.min():
    xmn += dx
xmx = np.round(x.max()/dx, 0)*dx
if xmx > x.max():
    xmx -= dx
ymn = np.round(y.min()/dy, 0)*dy
if ymn < y.min():
    ymn += dy
ymx = np.round(y.max()/dy, 0)*dy
if ymx > y.max():
    ymx -= dy

# Interpolate velocities to regular grid
X_unique = np.arange(xmn, xmx+dx/2, dx)
Y_unique = np.arange(ymn, ymx+dy/2, dy)
X, Y = np.meshgrid(X_unique, Y_unique)  # 2D grid for interpolation
interp = LinearNDInterpolator(list(zip(x, y)), v)
V = interp(X, Y)

# Extract depth of velocity isoline
# Go through every column and find first occurrence of v > desired velocity
# Then calculate the precise position through linear interpolation
if mode == 0:
    z_iso = []
    x_iso = []
    for i in range(len(X_unique)):
        n = np.where(V[:, i] > v_search)[0]
        if len(n) > 0:
            j = n[0]
            x_iso.append(X_unique[i])
            z_iso.append((Y_unique[j]-Y_unique[j-1])/(V[j, i]-V[j-1, i]) *
                         (v_search-V[j-1, i])+Y_unique[j-1])
        else:
            continue
# set limits of color bar scaling
    if v_min == v_max:
        v_min = np.floor(np.nanmin(V)/10.)*10.
        v_max = np.ceil(np.nanmax(V)/10.)*10.

# Plot velocity model and depth of isoline
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    pl = ax.imshow(V, cmap="rainbow", extent=[xmn, xmx, ymx, ymn],
                   vmin=v_min, vmax=v_max)
    ax.plot(x_iso, z_iso, "k")
    ax.set_xlabel("Distance[m]")
    ax.set_ylabel("Depth [m]")
    ax.set_title(f"{title} Isoline v={v_search:0.0f} m/s")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.2)
    cb = fig.colorbar(pl, cax=cax, label="velocity [m/s]")
    if scaling:
        ax.set_aspect("equal")
    fig.savefig(f"isoline_{v_search:0.0f}.png")
# Save isoline coordinates to file
    with open(f"isoline_{v_search:0.0f}.txt", "w") as fo:
        fo.write(f"x, z for isoline v={v_search:0.0f} m/s from file {file}\n")
        for i in range(len(x_iso)):
            fo.write(f"{x_iso[i]:0.1f}  {z_iso[i]:0.1f}\n")
# Extract vertical section
elif mode == 1:
    ncol = int((x_extract-X_unique[0])/dx)
    v_col = V[:, ncol]
# set limits of color bar scaling
    if v_min == v_max:
        v_min = np.floor(np.nanmin(V)/10.)*10.
        v_max = np.ceil(np.nanmax(V)/10.)*10.

# Plot velocity model and vertical section
    fig = plt.figure(figsize=(15, 6), tight_layout=True)
    # fig = plt.figure(figsize=(15, 6))
    # if scaling:
    #     fig.set_constrained_layout(True)
    gs = GridSpec(15, 13, figure=fig)
    ax_mod = fig.add_subplot(gs[:, :8])
    ax_sec = fig.add_subplot(gs[:, 10:])
    pl = ax_mod.imshow(V, cmap="rainbow", extent=[xmn, xmx, ymx, ymn],
                       vmin=v_min, vmax=v_max)
    ax_mod.plot([x_extract, x_extract], [ymx, ymn], "k", linewidth=3)
    ax_mod.set_xlabel("Distance[m]")
    ax_mod.set_ylabel("Depth [m]")
    ax_mod.set_title(f"{title} velocity model")
    if scaling:
        ax_mod.set_aspect("equal")
    if scaling:
        width = xmx - xmn
        height = ymx - ymn
        ax_mod.set_box_aspect(height / width)
    else:
        ax_mod.set_aspect("auto")        
    divider = make_axes_locatable(ax_mod)
    cax = divider.append_axes("right", size="2%", pad=0.2)
    cb = fig.colorbar(pl, cax=cax, label="velocity [m/s]")
    ax_sec.plot(v_col, Y_unique, linewidth=3)
    ax_sec.set_xlabel("Velocity [m/s]")
    ax_sec.set_ylabel("Depth [m]")
    ax_sec.set_title(f"Section at x={x_extract:0.1f} m")
    ax_sec.set_ylim((ymx, ymn))
    ax_sec.set_xlim((v_min, v_max))
    ax_sec.grid(axis="both")
    if scaling:
        ax_mod.set_ylim(ymx, ymn)
        ax_sec.set_ylim(ymx, ymn)

        # Keep both axes with identical physical y-size at draw time, even when
        # window size or layout engine changes positions.
        def _sec_locator(_ax, _renderer):
            pos_mod = ax_mod.get_position()
            pos_sec = ax_sec.get_position(original=True)
            return Bbox.from_bounds(
                pos_sec.x0, pos_mod.y0, pos_sec.width, pos_mod.height)

        ax_sec.set_axes_locator(_sec_locator)
    fig.savefig(f"vertical_section_{x_extract:0.0f}.png")
# Save y-coordinate and velocity to file
    with open(f"vertical_section_{x_extract:0.0f}.txt", "w") as fo:
        fo.write(f"z, velocities at x={x_extract:0.1f} m from file {file}\n")
        for i, v in enumerate(v_col):
            if np.isfinite(v):
                fo.write(f"{Y_unique[i]:0.1f}  {v:0.0f}\n")
# Extract horizontal section
else:
    nrow = int((y_extract-Y_unique[0])/dy)
    v_row = V[nrow, :]

# set limits of color bar scaling
    if v_min == v_max:
        v_min = np.floor(np.nanmin(V)/10.)*10.
        v_max = np.ceil(np.nanmax(V)/10.)*10.

# Plot velocity model and horizontal section
    fig = plt.figure(figsize=(10, 15), constrained_layout=True)
    gs = GridSpec(15, 100, figure=fig)
    ax_mod = fig.add_subplot(gs[:5, 3:90])
    ax_sec = fig.add_subplot(gs[7:, 3:90], sharex=ax_mod)
    ax_cb = fig.add_subplot(gs[:5, 92:95])
    pl = ax_mod.imshow(V, cmap="rainbow", extent=[xmn, xmx, ymx, ymn],
                       vmin=v_min, vmax=v_max)
    ax_mod.plot([xmn, xmx], [y_extract, y_extract], "k", linewidth=3)
    ax_mod.set_xlabel("Distance[m]")
    ax_mod.set_ylabel("Depth [m]")
    ax_mod.set_title(f"{title} velocity model")
    if scaling:
        width = xmx - xmn
        height = ymx - ymn
        ax_mod.set_box_aspect(height / width)
    else:
        ax_mod.set_aspect("auto")
    cb = fig.colorbar(pl, cax=ax_cb, label="velocity [m/s]")
    ax_sec.plot(X_unique, v_row, linewidth=3)
    ax_sec.set_xlabel("Distance [m]")
    ax_sec.set_ylabel("Velocity [m/s]")
    ax_sec.set_title(f"Section at z={y_extract:0.1f} m")
    ax_sec.grid(axis="both")
    ax_sec.set_ylim((v_min, v_max))
    if scaling:
        ax_mod.set_xlim(xmn, xmx)
        ax_sec.set_xlim(xmn, xmx)

        # Keep both axes with identical physical x-size at draw time, even
        # when window size or layout engine changes positions.
        def _sec_x_locator(_ax, _renderer):
            pos_mod = ax_mod.get_position()
            pos_sec = ax_sec.get_position(original=True)
            return Bbox.from_bounds(
                pos_mod.x0, pos_sec.y0, pos_mod.width, pos_sec.height)

        ax_sec.set_axes_locator(_sec_x_locator)
    fig.savefig(f"horizontal_section_{y_extract:0.0f}.png")
# Save x-coordinate and velocity to file
    with open(f"horizontal_section_{y_extract:0.0f}.txt", "w") as fo:
        fo.write(f"x, velocities at z={y_extract:0.1f} m from file {file}\n")
        for i, v in enumerate(v_row):
            if np.isfinite(v_row[i]):
                fo.write(f"{X_unique[i]:0.1f}  {v_row[i]:0.0f}\n")

# Keep figure open after computations; close window to finish script.
app = QtWidgets.QApplication.instance()
if app is not None:
    app.setQuitOnLastWindowClosed(True)
plt.show(block=True)
