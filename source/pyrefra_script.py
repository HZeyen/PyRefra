"""
Last modification Nov 25 2024

@author: Hermann Zeyen <hermann.zeyen@universite-paris-saclay.fr>
         Universite Paris-Saclay, France

Wrapper for pyrefra program
"""
import sys
import os
from pathlib import Path
from PyQt5 import QtWidgets
# sys.path.append("..")

from pyrefra import Pyrefra

if __name__ == "__main__":
    dir0 = os.getcwd()
    if isinstance(__package__, str):
        test2 = len(__package__) == 0
    else:
        test2 = __package__ is None

    def my_exception_hook(exctype, value, tracebk):
        """
        Test to capture CTLR-C, but does not work...

        Parameters
        ----------
        exctype : TYPE
            DESCRIPTION.
        value : TYPE
            DESCRIPTION.
        tracebk : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        print(exctype, value, tracebk)
        sys._excepthook(exctype, value, tracebk)
        sys.exit(1)

    if __name__ == "__main__" and test2:
        file = Path(__file__).resolve()
        parent, top = file.parent, file.parents[0]

        sys.path.append(str(top))
        try:
            sys.path.remove(str(parent))
# Already removed
        except ValueError:
            pass
        if top not in sys.path:
            sys.path.append(top)

        __package__ = "Pyrefra"

    try:
        app = QtWidgets.QApplication(sys.argv)
        sys._excepthook = sys.excepthook
        sys.excepthook = my_exception_hook
        main = Pyrefra.Main(str(top), dir0)
        main.window.showMaximized()
        sys.exit(app.exec_())
    except Exception as error:
        print(f'An unexpected exception occurred: {error}.')
        pass
