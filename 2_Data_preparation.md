# Data Preparation

The program is based on the files produced by DMT stations like Summit2 and Summit_X. These two systems have in common that (at least with our configuration), no coordinates are stored in the file / trace headers. Therefore, the program requires two **geometry files**:

**shots.geo**
	File contains four columns
    
		Shot point number (as recorded during acquisition)
        
		X_shot [m]
        
		Y_shot [m]
        
		Z_shot [m] (positive upwards)
        
In this file, every shot point number should have different coordinates associated. If due to recording errors two shots at the same coordinates were recorded with different shot point numbers, use file **file_corrections.dat** to correct shot point number for the corresponding file(s) (see below).

**receivers.geo**

This file contains the same kind of information as shots.geo, also in 4 columns. However, a fifth column may be added containing any letter describing the geophone component recorded in a trace (e.g., when working with 3-component geophones). You may select traces by this letter. Usually, letters are V or Z for vertical component, E or N for geographic directions, L or T for longitudinal or transverse components. However, any other letters are allowed. By default (no 5th column in a line), the letter “Z” is associated to a trace. So, you are not obliged to write the letter into all lines.

**If your system does not write the header words RECEIVER_STATION_NUMBER and SOURCE_STATION_NUMBER which are the ones corresponding to the first columns of files receivers.geo and shots.geo, the following definitions are made: RECEIVER_STATION_NUMBER is set to header keyword CHANNEL_NUMBER. SHOT_STATION_NUMBER is set to the number of your file (see below: NNNNN).** If this is not the good numbering, use file “file_corrections.dat” to change numbers. See below for explanations.
  
The program plots traces generally as function of offset which is calculated using the three coordinates. Therefore, in principle, absolute metric coordinates (e.g., UTM) may be given in both geometry files. However, we suggest strongly to convert coordinates to a local system and eliminate a possible trend, such that the line goes in X or Y direction. The program, internally, uses only the coordinate which has the largest extent (important for cdp points used, e.g, for distance gathers or export to SEGY format).

## Data file names
Files acquired with Summit equipment have by default the structure prefixNNNNN.ext. “prefix” is set by the user during acquisition. NNNNN is file number always with 5 ciphers filled to the left with zeros. “ext” may be “sg2” (Summit 2) or “seg2” (Summit X1). Also segy files are read with extension “sgy” or “segy”.

If the file names have a different structure, make sure that the file number is always just before the dot. The number of ciphers may be variable, not necessarily filled to the left with zeros.

Optionally, there may be two files that allow correction of certain acquisition errors:

## file_corrections.dat
This file is used to correct general errors for recorded files. You may only indicate files for which there are really problems. The file has six columns:

**irec, ipos, ir_start, ir_step, dt, interp**
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
This file is used to correct geophone problems occurring in all recorded shots. You may only indicate geophone positions where there is really a problem, no need to have a line for each geophone position It has two columns:
	**igeo, fact**

- **igeo**:	 number of geophone position (numbering starts at zero!)
- **fact**:	factor with which to multiply the data. This factor may be
	- 0 if the geophone recorded only noise
	- -1 if the polarity is wrong

If during interactive data treatment with program PyRefra.py amplitudes are changed or traces are muted, the program (when terminating normally) writes file **receivers_modify.dat**. This file may be renamed and used directly as receiver_corrections.dat.

PyRefra.py checks in the beginning if one or both of these files exist and if needed applies corrections in the data structure: ipos, ir_start and ir_step are used to change shot-point and receiver-point numbers in the headers, dt, fact and interp are directly applied to the corresponding trace data; if dt is not zero, the following procedure is applied: if dt is negative, `(-dt/sample-interval)` samples are taken away at the beginning of the traces and the same number of samples is added at the end with zero values. If dt is positive, `(dt/sample-interval)` samples are added at the beginning with zero values and the same number of samples is taken away at the end of the traces. In this way, the traces maintain their initial length in number of samples.


## PyRefra.config
Another optional file, **PyRefra.config**, may contain two text lines used for plot titles:
- title (used as title for most windows)
- direction contains the geographic direction of the beginning of the profile, may be one of
  		[S, SSW, SW, WSW, W, WNW, NW, NNW, N, NNE, NE, ENE, E, ESE, SE or SSE].
		The corresponding direction of the end of the profile is then determined automatically. So, if you enter e.g. SW for the beginning of the line, this means that the line runs in SW -> NE direction.
  
If this file does not exist at the beginning of a run, a dialogue box is opened asking for this information. The values are then stored in a newly created file PyRefra.config which will be used automatically the next time PyRefra is started.
