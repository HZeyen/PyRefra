import os
import sys
from PyQt5 import QtWidgets


class Files():
    """
    Contains methods linked to definition of files to be read
    """

    def __init__(self, dir0):
        self.dir = dir0
        self.numbers = []
        self.names = []
        self.file_dict = {}
        self.file_count = 0
        self.folder = ""
        self.file_ext = "seg2"
        self.file_type = "seg2"
        self.prefix = ""

    def check_names(self, files):
        """
        Check whether valid files have been chosen. If not, program finishes.
        Then check the list of chosen names and exclude all invalid ones
        (mainly this will be folder names)

        Parameters
        ----------
        files : list of strings
                Names of chosen items using
                QtWidgets.QFileDialog.getOpenFileNames from method get_files

        Returns
        -------
        files : list of strings
                Contains all valid file names

        """
        if len(files) == 0:
            print("No file chosen, program finishes")
            sys.exit("No file chosen")
        elif len(files[0]) == 0:
            print("\nNo file chosen, program finishes\n\nYou probably must "
                  + "close the Spyder console before restarting")
            sys.exit("No file chosen")
# Sort chosen file names
        fil = files[0]
        fil.sort()
# Check data format (SEG2 or SEGY)
        fname = fil[0]
        n = fname.rfind(".")
        self.file_ext = fname[n+1:]
        if self.file_ext in ("sg2", "seg2"):
            self.file_type = "seg2"
        elif self.file_ext in ("sgy", "segy"):
            self.file_type = "segy"
# If none of the above, stop program
        else:
            _ = QtWidgets.QMessageBox.critical(
                None, "Error",
                f"File type '{self.file_ext}' not recognized\n\n"
                + "Only sg2, seg2, sgy or segy allowed.\n\nProgram stops",
                QtWidgets.QMessageBox.Ok)
            raise NameError("file type error, should be sg2, seg2, sgy or "
                            + "segy\n")
# Loop over file names and exclude folders
        files = []
        for f in fil:
            if os.path.isdir(f):
                continue
            files.append(f)
        return files

    def get_number(self, file):
        """
        Determine coded number of recorded file
        Default file names from Summit2 instruments are PrefixNNNNN.sg2
        Default file names from Summit X1 instruments are PrefixNNNNN.seg2

        If not Summit files suppose that the file numbers are given just before
        the dot. If, by default, this is not the case, the file names should be
        changed before using PyRefra.py or corresponding files are ignored.

        Parameters
        ----------
        file : string
               File name to be analyzed

        Raises
        ------
        Exception
            Error if no number is found in file name.

        Returns
        -------
        num : int
              file number

        """
        n = file.rfind(".")
        num = -1
        for ipos in range(n):
            try:
                num = abs(int(file[ipos:n]))
                break
            except ValueError:
                continue
        if num < 0:
            answer = QtWidgets.QMessageBox.warning(
                None, "Warning",
                f"File {file} does not have standard numbering\n "
                + "Ignore this file or stop program and correct\n",
                QtWidgets.QMessageBox.Ignore | QtWidgets.QMessageBox.Close,
                QtWidgets.QMessageBox.Close)
            if answer == QtWidgets.QMessageBox.Close:
                raise ValueError("wrong file name, no number found")
        else:
            if ipos == 0:
                self.prefix = ""
            else:
                self.prefix = file[:ipos]
        return num

    def get_files(self):
        """
        Chose interactively files to be read, sort them alphabetically and
        extract folder name. Set this folder to working directory.
        Create dictionary file_dict containing file names and the numbers of
        the traces belonging to the corresponding file

        Returns
        -------
        None.

        """
# Open Explorer and show by default seg2 and sg2 files
        self.folder = None
        files = list(QtWidgets.QFileDialog.getOpenFileNames(
            None, "Select seismic data files", "",
            filter="seg2 (*.seg2 *.sg2) ;; segy (*.sgy *.segy) ;; all (*.*)"))
        files = self.check_names(files)

        self.folder = os.path.dirname(files[0])
        if len(self.folder) > 0:
            os.chdir(self.folder)
            print(f"folder: {self.folder}")
            print("Files read:")

        for nfil, f in enumerate(files):
            ff = os.path.basename(f)
# Default file names from Summit2 instruments are PrefixNNNNN.sg2
# Default file names from Summit X1 instruments are PrefixNNNNN.seg2
# Prefix is set by user during data acquisition
# If not Summit files suppose that the file numbers are given just before the
# dot if by default, this is not the case, the file names should be changed
# before using refraPy.py or corresponding files are ignored.
            num = self.get_number(ff)
            self.numbers.append(num)
            self.names.append(ff)
            self.file_dict[nfil] = {'name': ff}
            self.file_dict[nfil]['traces'] = []
        self.file_count = len(self.names)
