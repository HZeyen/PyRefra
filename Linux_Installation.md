The PyRefra has been installed and tested on a Ubuntu 22.04 plateform with the conda version 22.11.1


1. Open a terminal 
  + `conda update --all`
  + `conda config --add channels gimli --add channels conda-forge`
  + `conda create -n pg pygimli`
  + `conda activate pg`
  + `conda install obspy`
  + `conda install -c anaconda statsmodels`
  + `conda install -c anaconda spyder`
    
    
2. In **obspy**, you may find file **seg2.py** (for instance  in ~/anaconda3/envs/pg/lib/python3.8/site-packages/obspy/io/seg2/seg2.py)

There find line `if 'DELAY' in header['seg2']:` and comment out the whole "if" block. (a comment in python is the “#” character at the beginning of a line)

Also at the very end, find line containing `WARNING_HEADER` and comment it out.

This is not important for working of the program, but it avoids an annoyingly long warning message each time you start the program.

In **obspy < v.1.4**, there is a bug in reading seg2 files. Search in file **seg2.py** the line `if key == 'NOTE':` (line 337) in the line just above, change `value = ''` by `value = b''`

If under Linux, you get this error message: `AttributeError: 'numpy.int64' object has no attribute 'split'`, go to file **path_to_environment_site-packages/obspy/util/misc.py**, near line 217 and replace `except TypeError:` by simply `except:`

3. For nicer plots in **PyGimli**, modify numbering format in drawDataMatrix
The function is found in file **path_to_environment_site-packages/pygimli/viewer/mpl/dataview.py**

There, search lines starting with

`ax.set_xticklabels` and

`ax.set_yticklabels`

and change the rounding to 0 ciphers instead of the default 2 ciphers:

`ax.set_xticklabels(['{:g}'.format(round(xx[int(ti)], 0)) for ti in xt])`

4. In addition, it is interesting, though not absolutely necessary to add the following code to **pgimly** which allows summing up the number of rays having crossed each cell during all iterations and in this way plotting the final model avoiding the areas where no ray at all passed during any of the iteration steps:

File **path_to_environment_site-packages/pygimli/frameworks/inversion.py**

After line 522 (containing `startModel = self.startModel`) add:

`self.fop.cov_sum = np.zeros(len(self.fop.startModel()))`

After line 591 (starting with `self.modelHistory.append`) add:

`self.fop.cov_sum += self.fop.jacobian().transMult(np.ones(self.fop.jacobian().rows()))`



5. You can lauch spyder in the pg environment and run the example data-set
