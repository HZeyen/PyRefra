"""
Last modification Nov 25 2024

@author: Hermann Zeyen <hermann.zeyen@universite-paris-saclay.fr>
         Universite Paris-Saclay, France

Wrapper for pyrefra program. Allows starting the program interactively from
Spyder
"""

import sys
import os
from pathlib import Path
from PyQt5 import QtWidgets
# sys.path.append("..")

from pyrefra import Pyrefra

if __name__ == "__main__":
    dir0 = r"E:\Seg2Dat\Erreurs\UCT"

    os.chdir(dir0)

    def my_exception_hook(exctype, value, tracebk):
        """
        Not clear what it is good for.
        The problem is that usually no error messages are passed from
        QT.
        I found this in the internet, but it seems it does not work either...

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

    if isinstance(__package__, str):
        test2 = len(__package__) == 0
    else:
        test2 = __package__ is None
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
        print(f"An unexpected exception occurred: {error}.")
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)
        pass
