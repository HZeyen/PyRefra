# PyRefra
Software package for treatment and tomography of multishot refraction seismic data

## Installation
1. Mouse 

    It is strongly recommended to use a mouse. The program needs the central mouse button (or wheel). If you are working with a touchpad, configure it such that it simulates the central mouse button. The simplest is certainly to configure a three-finger-touch as central button. For this, under Windows press the WINDOWS-Key+I. On the screen that appears, click on “Périfériques” (the English version may be "Devices"). There, on the left-side panel, search “Pavé tactile” ("Touch Pad" in English?). Scroll down until you see the image with three fingers. Below this image, under “Appuis”, change to “Button central de la souris” ("Central mouse button"?).

2. Copy the following files all into the same folder of your choice: 

+ refraPy.py (the main program)
+ refraData.py
+ refraPlot.py
+ refraWindow.ui
3. Download and install Anaconda Individual Edition
4. Push the Windows key and search Anaconda3. There, choose (**RIGHT CLICK ON IT**) Anaconda Prompt and execute (if possible) as administrator. 
This opens a command window. Then type the following commands:
    + `conda update --all`
    + `conda config --add channels gimli --add channels conda-forge`
« channels » is the place where to find the source code, here « gimli » and “conda-forge”
    + `conda create -n pg pygimli pybert`
(this may take quite some time, don’t worry about error messages as long as conda is running). Here, “pg” means the new environment called “pg”. Pybert is not necessary for refraPy, but if you want also to use Orsay ERT inversion programs, better to install everything together.
        + `conda activate pg` (activation of the environment “pg”)
        + `conda install obspy`
        + `conda install statsmodels`
        + `conda update –all`  (just to be sure to have the latest versions)
        
5. Open Anaconda Navigator (Windows key -> Anadconda3 -> Anaconda Navigator)
    + In the upper tool bar, change “Applications on” from “base” or “anaconda” to “pg” (you will have to do this each time you open Anaconda!)
    + In the main window search for the icon “Spyder” and click on “Install” or, if it is already installed, on “Launch”.

6. Open Spyder
        + Open file refraPy.py (File -> open…)
        + In the Spyder tool bar click Run -> Configuration per file -> Execute in dedicated console

Then you may modify a few default settings in some packages for nicer plots or to avoid warnings:

The pg environment is installed at C:/Users/your_name/anaconda3/envs/pg/Lib (Linux: ~/anaconda3/envs/pg/lib/python3.8”). I will call this folder “ENV”

For nicer plots in **PyGimli**, modify numbering format in drawDataMatrix
The function is found in file **ENV/site-packages/pygimli/viewer/mpl/dataview.py**

There, search lines starting with

`ax.set_xticklabels` and

`ax.set_yticklabels`

and change the rounding to 0 ciphers instead of the default 2 ciphers:

`ax.set_xticklabels(['{:g}'.format(round(xx[int(ti)], 0)) for ti in xt])`

In addition, it is interesting, though not necessary to add the following code to **pgimly** which allows summing up the number of rais having crossed each cell during all iterations and in this way plotting the final model avoiding the areas where no rai at all passed during any of the iterations:

File **ENV/site-packages/pygimli/frameworks/inversion.py**

After line 522 (containing `startModel = self.startModel`) add:

`self.fop.cov_sum = np.zeros(len(self.fop.startModel()))`

After line 591 (starting with `self.modelHistory.append`) add:

`self.fop.cov_sum += self.fop.jacobian().transMult(np.ones(self.fop.jacobian().rows()))`

In **obspy**, you may find file **seg2.py** (usually in ENV/site-packages/obspy/io/seg2/seg2.py)

There find line `if 'DELAY' in header['seg2']:` and comment out the whole "if" block. (a comment in python is the “#” character at the beginning of a line)

Also at the very end, find line containing `WARNING_HEADER` and comment it out.

This is not important for working of the program, but it avoids an annoyingly long warning message each time you start the program.

If under Linux, you get this error message: `AttributeError: 'numpy.int64' object has no attribute 'split'`, go to file **ENV/site-packages/obspy/util/misc.py**, near line 217 and replace `except TypeError:` by simply `except:`
