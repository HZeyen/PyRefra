# Installation

1. Mouse 

    It is strongly recommended to use a mouse. The program needs the central mouse button (or wheel). If you are working with a touchpad, configure it such that it simulates the central mouse button. The simplest is certainly to configure a three-finger-touch as central button. For this, under Windows press the WINDOWS-Key+I. On the screen that appears, click on “Périfériques” (the English version may be "Devices"). There, on the left-side panel, search “Pavé tactile” ("Touch Pad" in English?). Scroll down until you see the image with three fingers. Below this image, under “Appuis”, change to “Button central de la souris” ("Central mouse button"?).

2. Download/clone the repository where the minimum required files are : 

+ PyRefra.py (the main program)
+ refraData.py
+ refraPlot.py
+ refraWindow.ui
+ PyRefra_Logo.png
3. Download and install Anaconda

4. Form here, the windows installation is presented, while [Linux Installations](Linux_Installation.md)  is provided as well.


5. Push the Windows key and search Anaconda3. There, choose (**RIGHT CLICK ON IT**) Anaconda Prompt and execute (if possible) as administrator. 
This opens a command window. Then type the following commands:
    1. `conda config --add channels gimli --add channels conda-forge`
    2. `conda install pygimli obspy statsmodels scikit-learn colorcet`
    3. `conda update --all`
   
   Sometimes conda blocks while installing. In this case, try to use 
   
    2. `conda install mamba`
    3. `mamba install pygimli obspy statsmodels scikit-learn colorcet`
    4. `mamba update --all`

   **If it is impossible to install pyGimli, you may continue without it, but the tomography option will not be available.**
        
6. Open Anaconda Navigator (Windows key -> Anadconda3 -> Anaconda Navigator)
    + In the upper tool bar, change “Applications on” from “base” or “anaconda” to “pg” (you will have to do this each time you open Anaconda!)
    + In the main window search for the icon “Spyder” and click on “Launch”.

7. Open Spyder
        + Open file PyRefra.py (File -> open…)
        + In the Spyder tool bar click Run -> Configuration per file -> Execute in an external system terminal

Then you may **modify a few things in some packages** for nicer plots or to avoid warnings. The corresponding modified files may be found on this github site within a folder with the same path structure as in the **pygimli** and **obspy** site-packages in Anaconda. There, for Windows, the pg environment is installed at C:/Users/your_name/anaconda3/Lib (Linux: ~/anaconda3/lib/python3.9). I will call this folder “ENV”. Site-packages are installed in “ENV/site-packages”

You may simply copy the folders pygimli and obspy into the conda site-packages folder. The other option, if you do not want to copy-paste the provided folders, is to modify the following files:

For nicer plots in **PyGimli**, modify numbering format in drawDataMatrix
The function is found in file **ENV/site-packages/pygimli/viewer/mpl/dataview.py**

There, search lines starting with

`ax.set_xticklabels` and

`ax.set_yticklabels`

and change the rounding to 0 ciphers instead of the default 2 ciphers:

`ax.set_xticklabels(['{:g}'.format(round(xx[int(ti)], 0)) for ti in xt])`

In addition, it is interesting, though not absolutely necessary to add the following code to **pgimly** which allows summing up the number of rays and ray lengths having crossed each cell during all iterations and in this way plotting the final model avoiding the areas where no ray at all passed during any of the iteration steps:

File **ENV/site-packages/pygimli/frameworks/inversion.py**

After line 546 (containing `startModel = self.startModel`) add:

`self.cov_sum = np.zeros(len(self.fop.startModel()))`

`self.cell_rays = []`

After line 628 (starting with `self.modelHistory.append`) add:

`self.cov_sum += self.fop.jacobian().transMult(np.ones(self.fop.jacobian().rows()))`

`self.cell_rays.append([])`

`for i in range(len(self.cov_sum)):`

`    n = len(np.where(np.array(self.fop.jacobian().col(i))!=0)[0])`
    
`    self.cell_rays[-1].append(n)`

In **obspy**, you may find file **seg2.py** (usually in ENV/site-packages/obspy/io/seg2/seg2.py)

There find line `if 'DELAY' in header['seg2']:` and comment out the whole "if" block. (a comment in python is the “#” character at the beginning of a line)

Also at the very end, find line containing `WARNING_HEADER` and comment it out.

This is not important for working of the program, but it avoids an annoyingly long warning message each time you start the program.

On the other hand, in the same file, we detected a problem which appears on acquisition systems using Anglophone date coding, i.e. mm/dd/yyyy instead of the expected dd/mm/yyyy. A work-around, only to get the program reading these data as well has been programmed in lines 183-188 with a try-except construction which gives the correct date if the coding is mm/dd/yyyy and if the day of month is larger than 12, but confounds day and month for smaller days:

`            try:`

`                day, month, year = int(date[0]), MONTHS[date[1].lower()], \`

`                    int(date[2])`

`            except:`

`                day, month, year = int(date[1]), MONTHS[date[0].lower()], \`

`                    int(date[2])`

In **obspy < v.1.4**, there is a bug in reading seg2 files. Search in file **seg2.py** the line `if key == 'NOTE':` (line 337) in the line just above, change `value = ''` by `value = b''`

If the program stops after having written to the screen “folder: …” followed by “files read” without a list of files, the most probable reason is that an old obspy version had been reinstalled and you must do the above corrections again!

This issue has been fixed within obspy (28/12/2022, bug fix #3178) but the corrected version is only installed for Python >=3.9.

If under Linux, you get this error message: `AttributeError: 'numpy.int64' object has no attribute 'split'`, go to file **ENV/site-packages/obspy/util/misc.py**, near line 217 and replace `except TypeError:` by simply `except:`
