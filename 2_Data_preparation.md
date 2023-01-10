# Data Preparation

The program is based on the files produced by DMT stations like Summit2 and Summit_X. These two systems have in common that (at least with our configuration), no coordinates are stored in the file / trace headers. Therefore, the program requires two **geometry files**:

**shots.geo**
	File contains four columns
    
		Shot point number (as recorded during acquisition)
        
		X_shot [m]
        
		Y_shot [m]
        
		Z_shot [m] (positive downwards)
        
In this file, every shot point number should have different coordinates associated. If due to recording errors two shots at the same coordinates were recorded with different shot point numbers, use file **file_corrections.dat** to correct shot point number for the corresponding file(s) (see below).

**receivers.geo**

	This file contains the same kind of information as shots.geo, also in 4 columns
  
The coordinates should be treated to eliminate a possible trend, such that the line goes in X or Y direction.

## Data file names
Files acquired with Summit equipment have by default the structure prefixNNNNN.ext. “prefix” is set by the user during acquisition. NNNNN is file number always with 5 ciphers filled to the left with zeros. “ext” may be “sg2” (Summit 2) or “seg2” (Summit X1).

If the file names have a different structure, make sure that the file number is always just before the dot. Number of ciphers may be variable, not necessarily filled to the left with zeros.

Optionally, there may be two files that allow correction of certain acquisition errors:

## file_corrections.dat
This file is used to correct general errors for recorded files. You may only indicate files for which there are really problems. The file has five columns:

**irec, ipos, ir_start, ir_step, dt**
- **irec**:	Number of record (file number)
DMT files have names of the following form:
prefixnnnnn.sg2 for Summit 2 or prefixnnnnn.seg2 for Summit_X
prefix may be user-defined, often “rec” or “shot_”
nnnnn is shot number (not shot-point number! The same shot point may be recorded in different files, if shots are repeated) in 5 ciphers, filled with zeros to the left.
- **ipos**:	Number of shot point position (counting starts with SP 1, not 0, i.e., natural numbering)
	If ipos is <= 0, actual shot point number is left unchanged
- **ir_start**:	number of first receiver position (natural numbering)
- **ir_step**:	step in numbering of receiver points (e.g., 2 if every second receiver point has been occupied)
	If ir_start and ir_step = 0, receiver station numbers remain unchanged
- **dt**:	trigger time correction [s]. If trigger is too late, i.e., the signal at 0m offset arrives before zero time, give a positive number, if signal arrives too late, give a negative number. If no geophone is placed at zero offset, you may use the air wave to check the trigger offset, perhaps after having applied a high-pass filter.
- **interp**	factor with which to increase the number of samples per trace (interpolation). I.e., if sampling interval is too large

## receiver_corrections.dat
This file is used to correct geophone problems occurring in all recorded shots. You may only indicate geophone positions where there is really a problem, no need to have a line for each geophone position It has three columns:
	**igeo, fact, interp**

- **igeo**:	 number of geophone position (numbering starts at zero!)
- **fact**:	factor with which to multiply the data. This factor may be
	- 0 if the geophone recorded only noise
	- -1 if the polarity is wrong

If during interactive data treatment with program refraPY.py amplitudes are changed or traces are muted, the program (when terminating normally) writes file **receivers_modify.dat**. This file may be renamed and used directly as receiver_corrections.dat.

refrapy.py checks in the beginning if one or both of these files exist and if needed applies corrections in the data structure: ipos, ir_start and ir_step are used to change shot-point and receiver-point numbers in the headers, dt, fact and interp are directly applied to the corresponding trace data; if dt is not zero, the following procedure is applied: if dt is negative, `(-dt/sample-interval)` samples are taken away at the beginning of the traces and the same number of samples is added at the end with zero values. If dt is positive, `(dt/sample-interval)` samples are added at the beginning with zero values and the same number of samples is taken away at the end of the traces. In this way, the traces maintain their initial length in number of samples.

