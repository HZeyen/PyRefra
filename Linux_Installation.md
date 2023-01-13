The PyRefra has been installed and tested on a Ubuntu 22.04 plateform with the conda version 22.11.1


1. Open a terminal 
  + `conda update --all`
  + `conda config --add channels gimli --add channels conda-forge`
  + `conda create -n pg pygimli`
  + `conda activate pg`
  + `conda install obspy`
  + `conda install -c anaconda statsmodels`
  + `conda install -c anaconda spyder`
    
    
2. In **obspy**, you may find file **seg2.py** (usually in /envs/sys/lib/pythonXX/site-packages/obspy/io/seg2/seg2.py)

There find line `if 'DELAY' in header['seg2']:` and comment out the whole "if" block. (a comment in python is the “#” character at the beginning of a line)

Also at the very end, find line containing `WARNING_HEADER` and comment it out.

This is not important for working of the program, but it avoids an annoyingly long warning message each time you start the program.

In **obspy < v.1.4**, there is a bug in reading seg2 files. Search in file **seg2.py** the line `if key == 'NOTE':` (line 337) in the line just above, change `value = ''` by `value = b''`

If under Linux, you get this error message: `AttributeError: 'numpy.int64' object has no attribute 'split'`, go to file **ENV/site-packages/obspy/util/misc.py**, near line 217 and replace `except TypeError:` by simply `except:`

3. You can lauch spyder in the pg environment and run the example data-set
