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

## Preparation

The program is based on the files produced by DMT stations like Summit2 and Summit_X. These two systems have in common that (at least with our configuration), no coordinates are stored in the file / trace headers. Therefore, the program requires two **geometry files**:

**shots.geo**
	File contains four columns
    
		Shot point number (as recorded during acquisition)
        
		X_shot [m]
        
		Y_shot [m]
        
		Z_shot [m] (positive downwards)
        
In this file, every shot point number should have different coordinates associated. If due to recording errors two shots at the same coordinates were recorded with different shot point numbers, use file “file_corrections.dat” to correct shot point number for the corresponding file(s) (see below).
receivers.geo
	File contains the same kind of information as shots.geo, also in 4 columns
The coordinates should be treated to eliminate a possible trend, such that the line goes in X or Y direction.
Data file names: Files acquired with Summit equipment have by default the structure prefixNNNNN.ext. “prefix” is set by the user during acquisition. NNNNN is file number always with 5 ciphers filled to the left with zeros. “ext” may be “sg2” (Summit 2) or “seg2” (Summit X1).
If the file names have a different structure, make sure that the file number is always just before the dot. Number of ciphers may be variable, not necessarily filled to the left with zeros.
Optionally, there may be two files that allow correction of certain acquisition errors:
file_corrections.dat
This file is used to correct general errors for recorded files. You may only indicate files for which there are really problems. The file has five columns:
irec, ipos, ir_start, ir_step, dt
irec:	Number of record (file number)
DMT files have names of the following form:
prefixnnnnn.sg2 for Summit 2 or prefixnnnnn.seg2 for Summit_X
prefix may be user-defined, often “rec” or “shot_”
nnnnn is shot number (not shot-point number! The same shot point may be recorded in different files, if shots are repeated) in 5 ciphers, filled with zeros to the left.
ipos:	Number of shot point position (counting starts with SP 1, not 0, i.e., natural numbering)
	If ipos is <= 0, actual shot point number is left unchanged
ir_start:	number of first receiver position (natural numbering)
ir_step:	step in numbering of receiver points (e.g., 2 if every second receiver point has been occupied)
	If ir_start and ir_step = 0, receiver station numbers remain unchanged
dt:	trigger time correction [s]. If trigger is too late, i.e., the signal at 0m offset arrives before zero time, give a positive number, if signal arrives too late, give a negative number. If no geophone is placed at zero offset, you may use the air wave to check the trigger offset, perhaps after having applied a high-pass filter.
interp	factor with which to increase the number of samples per trace (interpolation). I.e., if sampling interval is too large

receiver_corrections.dat
This file is used to correct geophone problems occurring in all recorded shots. You may only indicate geophone positions where there is really a problem, no need to have a line for each geophone position It has three columns:
	igeo, fact, interp

igeo:	 number of geophone position (numbering starts at zero!)
fact:	factor with which to multiply the data. This factor may be
	0 if the geophone recorded only noise
	-1 if the polarity is wrong
interp	factor with which to increase the number of samples per trace (interpolation). I.e., if sampling interval is too large
If during interactive data treatment with program refraPY.py amplitudes are changed or traces are muted, the program (when terminating normally) writes file “receivers_modify.dat”. This file may be renamed and used directly as receiver_corrections.dat.
Seg2_read.py checks in the beginning if one or both of these files exist and if needed applies corrections in the data structure: ipos, ir_start and ir_step are used to change shot-point and receiver-point numbers in the headers, fact and interp are directly applied to the corresponding trace data; if dt is not zero, the following procedure is applied: if dt is negative, (-dt/sample-interval) samples are taken away at the beginning of the traces and the same number of samples is added at the end with zero values. If dt is positive, (dt/sample-interval) samples are added at the beginning with zero values and the same number of samples is taken away at the end of the traces; In this way, the traces maintain their initial length in number of samples.

