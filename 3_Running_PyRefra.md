# Running PyRefra

To start **refrapy.py**, the best is to open it in Spyder and click on the green arrow or press F5. **If you configured Spyder as indicated in “Installation.5”, and you restart refraPy, you must first close the console of the earlier run. If not, a possible change of the folder will not be taken into account.**

Before running the program, $\textcolor{red}{\text{you must define the path to the python files}}$ and may want to set manually the path to your data folder. For this go to **approximately line 90**. There, you will see a line starting with `sys_path = r”…”`. Replace the existing path by the path of the python files. In the following line, starting with `self.dir0 = r”…”`, replace the existing path by the path to the data files. This is not necessary, since you may choose the path interactively, but if you work for a longer time on the same data set, it avoids searching on each program start your data on the disk.

The program will ask for the files to be opened. It proposes initially files with ending **seg2**. If none is shown in the dialog box, probably you are dealing with Summit2 or SEGY files, so you may change the default ending to **sg2** or **sgy** (or even show all files). Usually, you will want to open all available files. Click on one of them and chose them with **CTRL-A**. If not, you may choose them as usual in the explorer, using SHFT or CTRL keyboard keys. No problem, if also folders are chosen with CTRL-A, the program filters them out. If you know that there are bad files, you may copy them into another folder or change their name-endings to something else than “.seg2”, “.sg2” or .”sgy” in which case the default option will not find these files and they are not opened. File numbering does not need to be contiguous. $\textcolor{red}{\text{For information}}$: if segy files are read, internally, the program copies the important segy header entries to seg2 headers. If the program has problems, it might be that not all necessary header entries exist or have been copied over.

In the beginning, the program shows a warning:

`C:\Users\Hermann\anaconda3\envs\pg\lib\site-packages\pkg_resources\__init__.py:123: PkgResourcesDeprecationWarning: -PKG-VERSION is an invalid version and will not be supported in a future release`

This warning comes from “import obspy”. Hope, obspy developers will do something early enough…

When opened, the GUI shows in the main window the data of the shot gather corresponding to the first shot point. A smaller window to the right shows all available shot point numbers. Clicking on one of the names plots the seismograms of the corresponding shot point in the main window. Beneath the main window, a help text is shown. To be visible, it may be necessary that the Windows task bar be automatically hidden.

![Main_window_screen_copy](https://github.com/HZeyen/PyRefra/blob/main/Screen_shot_refrapy.png)

