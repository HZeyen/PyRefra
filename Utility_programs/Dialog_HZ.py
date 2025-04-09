# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 12:41:31 2024

@author: Hermann
"""

from PyQt5 import QtWidgets, QtCore

class Dialog(QtWidgets.QWidget):
    """
    Created on Tue Sep 13 12:50:26 2022
    @author: Hermann
    General purpose dialog box

    Input:
        parent (Class object): Certain values are needed from the parent class
                               usually, the call will thus be:
                               Dialog(self,labels,types,values,title)
        labels (list, string): labels explaining the different fields
                               For a series of radiobuttons, labels[i] is itself
                               a list
        types (list,string): In the order of "labels", indicates the type
                            of field the label belongs to. If labels[i] is itself
                            a list, this is considered as one single label.
                            May be
                            'l' for label (no dialog entry, only information)
                            'c' for checkbox
                            'e' for LineEdit (to enter values)
                            'r' for radiobutton
                            'b' for combo box (drop down menu)
                            may be capital or not
        values (list, string, float or int): initial values for LineEdit fields,
                            number of radiobutton activated by default.
                            for labels, it may be "b" (bold), "i" (italic), or
                            anything else (usually None) for standard text
                            Ignored for checkbox.
                            Optional, default value: None
        title (string): Title of the dialog window
                            Optional, default value: "Title"
    """
# Create window
    def __init__(self, labels, types, values = None, title="Title"):
        import numpy as np
        super(Dialog, self).__init__()
        nlab = len(labels)
        self.labels = labels
        self.types = types
        self.ck_order = np.zeros(nlab, dtype = int)
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
        for i in range(nlab):
            self.ckState.append(False)
        self.YesBtn = QtWidgets.QPushButton('Ok',self)
        self.YesBtn.move(10,20*(nlab+3))
        self.CancelBtn = QtWidgets.QPushButton('Cancel',self)
        self.CancelBtn.move(150,20*(nlab+3))
        self.mainLayout = QtWidgets.QGridLayout()
        self.setLayout(self.mainLayout)
        il_add = 0
        ilin = 0
        for i in range(len(labels)):
            il = ilin + il_add
            if self.types[i].lower() == 'l':
                if values[i]:
                    if values[i].lower() == 'b':
                        self.dlabels.append(QtWidgets.QLabel("<b>"+labels[i]+"</b>"))
                    elif values[i].lower() == 'i':
                        self.dlabels.append(QtWidgets.QLabel("<i>"+labels[i]+"</i>"))
                    else:
                        self.dlabels.append(QtWidgets.QLabel(labels[i]))
                else:
                    self.dlabels.append(QtWidgets.QLabel(labels[i]))
                self.dlines.append(None)
                self.ckb.append(None)
                self.rbtn.append(None)
                self.btngroup.append(None)
                self.combo.append(None)
                self.mainLayout.addWidget(self.dlabels[ilin], il, 0, 1, 2)
                ilin += 1
            elif self.types[i].lower() == 'e':
                self.dlabels.append(QtWidgets.QLabel(labels[i]))
                self.dlines.append(QtWidgets.QLineEdit())
                self.ckb.append(None)
                self.rbtn.append(None)
                self.btngroup.append(None)
                self.combo.append(None)
                self.mainLayout.addWidget(self.dlabels[ilin], il, 0, 1, 1)
                self.mainLayout.addWidget(self.dlines[ilin], il, 1, 1, 1)
                try:
                    s = str(values[i])
                    self.dlines[-1].setText(s)
                except:
                    pass
                ilin += 1
            elif self.types[i].lower() == 'r':
                self.ckb.append(None)
                self.combo.append(None)
                self.rbtn.append([])
                self.btngroup.append(QtWidgets.QButtonGroup())
                rck = int(values[i])-1
                if rck<0 or rck>=len(labels[i]):
                    rck = 0
                for ir,l in enumerate(labels[i]):
                    self.dlabels.append(None)
                    self.dlines.append(None)
                    self.rbtn[i].append(QtWidgets.QRadioButton(l))
                    self.btngroup[-1].addButton(self.rbtn[i][-1])
                    self.mainLayout.addWidget(self.rbtn[i][-1], il, 0, 1, 2)
                    if ir == rck:
                        self.rbtn[i][-1].setChecked(True)
                    else:
                        self.rbtn[i][-1].setChecked(False)
                    il += 1
                    ilin += 1
            elif self.types[i].lower() == 'c':
                self.dlabels.append(None)
                self.dlines.append(None)
                self.rbtn.append(None)
                self.btngroup.append(None)
                self.combo.append(None)
                self.ckb.append(QtWidgets.QCheckBox(self))
                self.ckb[i].setText(self.labels[i])
                self.mainLayout.addWidget(self.ckb[i], il, 0, 1, 2)
                self.ckb[i].stateChanged.connect(self.checked)
                ilin += 1
            elif self.types[i].lower() == 'b':
                self.dlabels.append(None)
                self.dlines.append(None)
                self.ckb.append(None)
                self.rbtn.append(None)
                self.btngroup.append(None)
                self.combo.append(QtWidgets.QComboBox())
                for il,l in enumerate(labels[i]):
                    self.combo[i].addItem(l)
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
# If a check box is checked, search the one which has been checked
        if checked == Qt.Checked:
# Do necessary changes
            for i,ck in enumerate(self.ckb):
# If self.ckb has None value, the corresponding entry is not a checkbox
                if not ck:
                    continue
# If self.ckb.checkState is checked after click, set ckState to True and do
#    changes
                if ck.checkState()==Qt.Checked:
                    self.ckState[i] = True
# If checkbox nr i was not checked, increase the number of checked boxes
#    (n_checked) and store the order of checkin in self.ck_order
                    if self.ck_order[i] == 0:
                        self.n_checked += 1
                        self.ck_order[i] = self.n_checked
# if checkboxis checked, change the label, indicating the order in which it was
#    checked.
                        self.ckb[i].setText(f"{self.labels[i]} ("+\
                                            f"{self.ck_order[i]})")
                        break
# If self.ckb.checkState is still unchecked, set ckState to Falsee
                else:
                    self.ckState[i] = False
# If click has unchecked a checkbox, do necessary changes
        else:
            for i,ck in enumerate(self.ckb):
# If self.ckb has None value, the corresponding entry is not a checkbox
                if not ck:
                    continue
# If self.ckb.checkState is still checked, set ckState to True
                if ck.checkState()==Qt.Checked:
                    self.ckState[i] = True
# If checkbox is no longer checked but it was (self?ckState), the unchecked box
#    is found
                else:
                    if self.ckState[i] == True:
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
                                self.ckb[j].setText(f"{self.labels[j]} ("+\
                                                   f"{self.ck_order[j]})")
                        self.n_checked -= 1
                        break
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
              'r': to define a series of radiobuttons (the corresponding list of
                   labels is considered as one single label for the numbering of
                   labels, types and values)
              "b": to define a combobox (dropdown menu)
      values : list of texts, numbers or Nones (length as labels)
          Initial values given to editable fields. Should be None for check boxes
              should be the number (natural numbering, starting at 1, not at 0)
              of the radiobutton to be activated by default.
              For labels, it may be "b" (bold text), "i" (italic) or anything
              else, including None for standard text
      title : str, default: "Title"
          Title of the dialog box.

      Returns
      -------
      results : list of bool
          Response of each data entrance field. Should transformed to the
          needed data format (int, float...) in the calling function
  		For radiobuttons, the returned value indicates the number of the
  		active button (counting starts at 0). For checkboxes, the returned
  		value is -1 if the box was not checked and the position at which
  		the box was checked (starting at 0) if it has been checked
      Dbutton: bool
          If True, "Apply" button has been pressed to finish dialog, if False
          "Cancel" button has been pressed.
          
    """
    from PyQt5.QtWidgets import QApplication
    from PyQt5 import QtCore
    import sys
    _ = QApplication(sys.argv)
    D = Dialog(labels, types, values, title)
    D.Dfinish = False
    while (D.Dfinish != True):
        QtCore.QCoreApplication.processEvents()

    results = []
    iline = 0
    l = len(values)
    if l > 0:
        for it,t in enumerate(types):
            if t.lower() == "e":
                results.append(D.dlines[iline].text())
                iline += 1
            elif t.lower() == "r":
                results.append(None)
                for i in range(len(labels[it])):
                    iline += 1
                    if D.rbtn[it][i].isChecked():
                        results[-1] = i
#                            break
            elif t.lower() == 'c':
                results.append(D.ck_order[it]-1)
                iline += 1
            elif t.lower() == "b":
                results.append(D.combo[it].currentIndex())
                iline += 1
            else:
                results.append(None)
                iline += 1
    return results, D.Dbutton

