# Running PyRefra

To start **refrapy.py**, the best is to open it in Spyder and click on the green arrow or press F5. **If you configured Spyder as indicated in “Installation.5”, and you restart refraPy, you must first close the console of the earlier run. If not, a possible change of the folder will not be taken into account.**

Before running the program, $\textcolor{red}{\text{you must define the path to the python files}}$ and may want to set manually the path to your data folder. For this go to **approximately line 90**. There, you will see a line starting with `sys_path = r”…”`. Replace the existing path by the path of the python files. In the following line, starting with `self.dir0 = r”…”`, replace the existing path by the path to the data files. This is not necessary, since you may choose the path interactively, but if you work for a longer time on the same data set, it avoids searching on each program start your data on the disk.

The program will ask for the files to be opened. It proposes initially files with ending **seg2**. If none is shown in the dialog box, probably you are dealing with Summit2 or SEGY files, so you may change the default ending to **sg2** or **sgy** (or even show all files). Usually, you will want to open all available files. Click on one of them and chose them with **CTRL-A**. If not, you may choose them as usual in the explorer, using SHFT or CTRL keyboard keys. No problem, if also folders are chosen with CTRL-A, the program filters them out. If you know that there are bad files, you may copy them into another folder or change their name-endings to something else than “.seg2”, “.sg2” or .”sgy” in which case the default option will not find these files and they are not opened. File numbering does not need to be contiguous. $\textcolor{red}{\text{For information}}$: if segy files are read, internally, the program copies the important segy header entries to seg2 headers. If the program has problems, it might be that not all necessary header entries exist or have been copied over.

In the beginning, the program shows a warning:

`C:\Users\Hermann\anaconda3\envs\pg\lib\site-packages\pkg_resources\__init__.py:123: PkgResourcesDeprecationWarning: -PKG-VERSION is an invalid version and will not be supported in a future release`

This warning comes from “import obspy”. Hope, obspy developers will do something early enough…

When opened, the GUI shows in the $\textcolor{blue}{\text{main window}}$ the data of the shot gather corresponding to the first shot point. A smaller window to the right shows all available $\textcolor{red}{\text{shot point numbers}}$. Clicking on one of the names plots the seismograms of the corresponding shot point in the $\textcolor{blue}{\text{main window}}$. Beneath the main window, a $\textcolor{green}{\text{help text}}$ is shown. To be visible, it may be necessary that the Windows task bar be automatically hidden. Above the main window, a $\textcolor{violet}{\text{menu bar}}$ is placed.

![Main_window_screen_copy](https://github.com/HZeyen/PyRefra/blob/main/Screen_shot_refrapy.png)

The first thing to do will be to increase the GUI to full screen size. If the window is too small, certain graphs containing multiple partial windows (e.g., tomography results) will not show up correctly, even in the stored png files (overlapping axis texts of different sub-windows).

## $\textcolor{red}{\text{Actions without Menu:}}$

**Numpad or standard + or -**

Increase or decrease amplitudes of traces

**Keyboard arrows right or left**
If you have made a zoom in X (showing only part of the traces), the arrow right will shift the view to the traces further to the right (larger X coordinates), the arrow to the left will show traces with smaller X coordinates. The size of the window in terms of X-coordinates will always be maintained. E.g., your zoom shows traces 1 to 50 out of 120 traces. The Right-Key will show traces 51 to 100, a further click will show traces 71 to 120 (still 50 traces are shown, but since the last trace is n° 120, the window starts at 71). A Left-Key click will now show traces 21 to 70 and a further Left-Key click will show again traces 1-50.

## The Menus

### **$\textcolor{red}{\text{File Menu}}$**

**$\textcolor{violet}{\text{Choose Data File}}$** (keyboard shortcut: **CTRL-O**): 

 Usually, you will not need this menu point. It allows you to choose new data sets.

**$\textcolor{violet}{\text{Save SEGY}}$** (keyboard shortcut: **Y**): 

 Save data in SEGY format (e.g., to continue with Seismic Unix)
 A dialog box opens where you may choose a series of options:
   - Data to be saved: 
     - only data visible on the screen (if zoomed)
     - all data, including possible data recorded before the trigger
     - all data recorded after the trigger
   - Data of all shots or only of the active shot
   - Save all data in one single file (this file will be called prefix00000.sgy) or each shot in another file

**$\textcolor{violet}{\text{Save SU}}$** (keyboard shortcut: **U**): 

Save data in SU format (e.g., to continue with Seismic Unix). The same dialog box as for SEGY storing appears. It seems that in SU format, several header words are not correctly filled in. Therefore, in my experience, it is better to save the data in SEGY format and use the following command n SU (Linux) in order to erase the file header and pass to SU format:
`segyread tape=file_name verbose=1 | segyclean > data.su`

**$\textcolor{violet}{\text{Save SEG2}}$** (keyboard shortcut: **2**): 
Save data in SEG2 format (e.g., if you have done quite a few corrections via file_corrections.dat or receiver_corrections.dat and/or a number of files were eliminated, and you want a new list with continuous file numbering). You may save just the gather on the screen or all available gathers (into a single file or every gather into its own file). You may write file-gathers or shot-gathers. If only the gather on the screen is to be saved, also only the corresponding type of gather is written, where receiver gathers and distance gathers cannot be stored, and a warning message appears. File gathers are written into the interactively chosen folder as “file_nnnnn.seg2” and shot-gathers as “shot_nnnnn.seg2”, nnnnn being the renumbered file number (starting at 1) or the shot number.

**$\textcolor{violet}{\text{Save FWI}}$** (no keyboard shortcut): 
Save data only in binary format without any header (4 bytes)

**$\textcolor{violet}{\text{Save ASCII}}$** (no keyboard shortcut): 
Save data in ASCII format
- Line 1: t0, dt, ntrace, nsamp, maximum value, factor, x_source [m], y_source [m], z_source [m]
-- t0: time of first sample [s]; dt: sampling interval [s]; ntrace: number of traces; nsamp: number of samples per trace
- Line 2: ntrace columns containing the value of each trace corresponding to “maximum value” of line 1
The amplitudes are calculated as stored_value/(maximum_value*factor)*corresponding_value
“stored_value”: Integer values stored in the rest of the file; “corresponding_value” value in this Line 2 for each trace
- Lines 3 to 5: X_receivers, Y_receivers, Z_receivers (each line ntrace values)
- Following nsamp lines: Normalized values of all traces (one line per time interval)

**$\textcolor{violet}{\text{Save Plot}}$** (keyboard shortcut: **P**): 
Save actual screen into png file. File names depend on what is actually shown

**$\textcolor{violet}{\text{Quit}}$** (keyboard shortcut: **CTRL-Q**): 
Finish program with cleaning up. A dialog window opens to confirm termination.

### **$\textcolor{red}{\text{Plot Menu}}$**

**$\textcolor{violet}{\text{Original data this screen}}$** (keyboard shortcut: **O**): 
Plot the original, not treated data shown on the actual screen (e.g., to undo filtering). All other data are not reset. A potential zoom is maintained.

**$\textcolor{violet}{\text{Original data all}}$** (keyboard shortcut: **SHFT-O**): 
Plot the original, not treated data of all files (e.g., to undo airwave muting or velocity filter applied to all files). A potential zoom is maintained.

**$\textcolor{violet}{\text{Shot gather}}$** (keyboard shortcut: **S**): 
Plot shot gathers (in contrast to file gather), i.e., all traces from all acquired files that have been recorded from the same shot point position are plotted into one record section. If the same combination receiver point – shot point has been recorded several times, the traces are plotted overlapping. The list of plots on the “Available shots” window is updated and the gather corresponding to the first shot is shown in the main window.

**$\textcolor{violet}{\text{Plot data of one file}}$** (keyboard shortcut: **F**): 
Plot file gathers (in contrast to shot gather), i.e., the traces from the acquired files (one file = one shot) are plotted into one record section. The list of plots on the “Available shots” window is updated and the gather corresponding to the first file is shown in the main window.

**$\textcolor{violet}{\text{Receiver gather}}$** (keyboard shortcut: **R**): 
Plot receiver gathers, i.e., all traces from all acquired files that have been recorded at the same receiver position are plotted into one record section. If the same combination receiver point – shot point has been recorded several times, the traces are plotted overlapping. The list of plots on the “Available shots” window is updated and the gather corresponding to the first receiver is shown in the main window.

**$\textcolor{violet}{\text{Distance gather}}$** (keyboard shortcut: **D**): 
Plot distance gathers, i.e., all traces from all acquired files that have a common absolute offset are plotted into one record section. If the same offsets have been recorded several times, the traces are plotted with a slight lateral shift. The traces are plotted at the mid-point position between shot and receiver. The list of plots on the “Available shots” window is updated and the gather corresponding to the smallest offset is shown in the main window. This window may serve for quality control of the picks. Especially, in the case when receiver and shot point are switched both traces are plotted at the same position and should have the same arrival time.

**$\textcolor{red}{\text{Zoom options}}$**:
Most zoom options are maintained from one plotting window to the others

**$\textcolor{violet}{\text{Zoom}}$** (keyboard shortcut: **Z**): 
Trace a rectangle with left mouse key maintained pressed to select data to be shown (as well in time as in distance). When the mouse button is released, the zoom is accepted. The time limits of the zoom will be maintained from one plot to the next, however not the x limits. When calling a new plot, again all traces are plotted.

**$\textcolor{violet}{\text{Zoom out}}$** (keyboard shortcut: **CTRL-Z**): 
If several zoom levels have been defined using “Zoom”, return to the next smaller zoom.

**$\textcolor{violet}{\text{Zoom in}}$** (keyboard shortcut: **SHIFT-Z**): 
If several zoom levels have been defined using “Zoom” and you made a Zoom-out (or Zoom-initial), return to the next tighter zoom.

**$\textcolor{violet}{\text{Zoom initial}}$** (keyboard shortcut: **ALT-Z**): 
Go back to the original plotting limits, i.e., show all data.

**$\textcolor{red}{\text{Gain options}}$**:
Gain options are maintained from one plotting window to the others.

**$\textcolor{violet}{\text{Trace normalize}}$** (no keyboard shortcut): 
Plot the traces normalized to the maximum value of each trace (default setting).

**$\textcolor{violet}{\text{Time gain}}$** (no keyboard shortcut): 
Apply a time gain (traces are in addition trace normalized). A dialog window opens, asking for the exponent of the time gain, i.e., amplitudes are multiplied by (t**time_gain).

**$\textcolor{violet}{\text{Distance gain}}$** (no keyboard shortcut): 
Apply a distance gain (traces are not trace normalized). A dialog window opens, asking for the exponent of the distance gain, i.e., amplitudes are multiplied by (x**distance_gain).

**$\textcolor{violet}{\text{AGC}}$** (no keyboard shortcut): 
Apply an Automatic Gain Control (traces are in addition trace normalized). A dialog window opens, asking for the AGC window length.

### **$\textcolor{red}{\text{Utilities Menu}}$**

**$\textcolor{violet}{\text{P Model}}$** (keyboard shortcut: **CTRL-P**): 
Allows to construct a 1D velocity model by tracing a series of straight lines across the record section.

To trace a line, press left mouse button and pull the line. Releasing the mouse button is accepting the line.

The first straight line is considered to be the direct wave and is forced to pass through the origin. For all further lines, the position where the mouse is clicked is one point of the line, the position where the mouse is released is the second point.

When releasing the mouse, a dialog window opens which allows you to accept the line and trace the next one, to accept it and finish the model or to redo the line. In addition, when clicking on “Show details”, the velocity, and intercept times are given as well as the depth to the layer limit (refractor) and the thickness of the layer above the refractor.
The lines plotted will disappear when a new plot is called (right window or different gather). A new call to P_Model creates a new model, starting with the direct wave.
When the model is finished, it is written to a file in the data folder with the following name structure: prefix_number_data_time.1Dmod where prefix may be “shot_point”, “receiver_point” or “file” depending on the type of gather used (it is supposed that the routine is not used when a distance gather is shown.

**$\textcolor{violet}{\text{Tomography}}$** (keyboard shortcut: **T**): 
This tool is only activated if measured travel times exist. It uses **Pygimli** to invert travel times into a tomography model.
On call, a dialog box is opened asking for:
- Maximum depth of the model [m]: The default value is calculated as 1/3 of the maximum available offset
- Initial smoothing parameter: The bigger the value, the smoother will be the model. This will be a parameter to play with, doing different inversions with different smoothing
- Smoothing reduction per iteration. If not zero, the smoothing parameter will be multiplied by this factor after each iteration
- Smoothing in z direction: This value may be between 0 and 1. When 0, the different layers of the tomography will be independent from each other, if 1, the same smoothing factor will be used in vertical and horizontal direction.
- Maximum iterations: The number of iterations may be limited with this value, e.g., for testing purposes. In general, the iterations are stopped if the data fit does not become better anymore, i.e., if the chi² value does not decrease by more than 1% of its actual value or if chi² becomes smaller than 1.
- Initial and final velocities: Ray tracing needs an initial model with a vertical gradient. This initial model is defined by these two values. **TODO**: allow for a refined starting model.
- Minimum and maximum allowed velocities: Pygimli allows to limit search space of model velocities. Give here the desired extreme velocities. **TODO**: Understand working of Pygimli: Even if very large upper velocity limit is given, pygimli limits the velocity model to a range of values much narrower than given. Only if the parameter “limits=[vmin,vmax]” is not at all used, velocities seem to be completely free. Therefore, if both given velocities are zero, program makes a call to mgr.invert without using key word limits
- velocity color scale min and max: If you do inversions for different profiles, it may be interesting to have the same color scale for all. In this case, you may give the values of the minimum and maximum velocities appearing in the color bar.
- Plot title: The default title is the name of the data folder. You may change it to any text. Appears only at main title of the inversion results plot.

At the end of the inversion, the resulting model is shown in the upper 60% of the screen. Below, smaller sub-windows show the initial model, the rays of the final model with the “coverage” (i.e., a measure of resolution), the measured travel times as function of shot point position and receiver point position, the misfits in a similar plot calculated as “calculated times minus measured times” and finally, a plot showing average misfits for the different shot points and the different receiver points. This last plot allows you to verify if there are problems with certain shots or receivers.
This plot is stored in a png file, together with other inversion results such as velocity model, misfits, rays of final model and chi²-evolution in a folder created automatically by Pygimli the name of which has the following structure: ./date-time/TravelTimeManager (where the “.” Is the data folder).
To leave this plot, click on an entrance in the “Available plots” window. To show the plot again, you should open the saved png outside of the program or you must do the inversion again.

**$\textcolor{violet}{\text{Tau-P}}$** (keyboard shortcut: **CTRL-T**): 
Tool calculates and shows Tau-P analysis. For the moment, only the plot is shown. If you want to save it, use $\textcolor{violet}{\text{File -> Save plot}}$

**$\textcolor{violet}{\text{False colour}}$** (keyboard shortcut: **SHFT-F**): 
Possibility to plot 1 to 3 different indicators that may help for picking. This tool is for the moment only meant for display, to analyse data visually. Possible indicators are: Instant frequency, envelope, 2nd derivative of envelope, 2nd derivative of data, Sta-Lta transform, Akaike Information Criterium (AIC), max-min relation and autocorrelation. Max-min relation plots the relative amplitude difference between a maximum and its following minimum on the one hand and the same maximum and its preceding minimum on the other hand. Autocorrelation is usually plotted on its own since the time scale has no relation with the other options. On top of the autocorrelation, two traces are plotted: white is the average autocorrelation, yellow the average trace obtained from averaging the frequency spectra. A dialog box will open where you may check the desired indicators (1, 2 or 3 out of the 8 possibilities). The first indicator you click on will determine the red channel, the second one the green channel and the third one the blue channel. This implies also that if only one indicator is chosen, the plot will be between black and red. If 2 are chosen, colour may change between black, red, green, and yellow (no blue). If picks have been measured for the analysed gather, they are plotted on top of the false colour plot. In addition, a magenta line is plotted indication evolution of maxima of the combined indicators (maximum brightness, Pythagoras of color values), as well as a white line which represents an optimum fit of the positions of the maxima through two straight lines. The maxima are possible picks, the straight line indicates a possible two-layer velocity model. These lines are though not always meaningful.

**$\textcolor{violet}{\text{Change trace sign}}$** (keyboard shortcut: **I**): 
You may choose one or more traces which have wrong polarity (geophone or receiver electronics may have been reverse-cabled). Left click on one trace after the other. The last trace to be chosen should be marked a small red star will be plotted near the bottom of the trace plot. When all traces have been marked, press right mouse key.
When the program is finished using $\textcolor{violet}{\text{File -> Quit}}$, a file “receivers modify.dat” is written indicating the muted traces, which may eventually be renamed to “receiver_corrections.dat” to make these changes durable.

**$\textcolor{violet}{\text{Change sign}}$** (keyboard shortcut: **CTRL-I**): 
Depending on the polarity of the acquisition system, the general polarity may be reversed. In order to show first arrivals in positive direction (filled black), you may use this tool that changes the sign of all recorded traces, not only those visible in the actual window.

**$\textcolor{violet}{\text{Change colors tomo}}$** (keyboard shortcut: **C**): 
This option is only activated when a tomography model has been calculated and is visible on the screen. Allows changing color scale of model and depth extent shown.

**$\textcolor{violet}{\text{Wave animation}}$** (keyboard shortcut: **W**): 
Plots animation of wave evolution in time along the geophone line actually on the screen. Animation goes from first to last sample of actual time zoom. Data are smoothed with 3-point running average and amplitudes are cut at 99% quantile.

**$\textcolor{violet}{\text{Attenuation}}$** (keyboard shortcut: **Q**): 
Function searches for each trace the maximum of the envelopes of data plotted on the screen (you may use muting functions to focus surgically on certain phases). Amplitudes are multiplied by the absolute offset to counteract geometric spreading. Then, for each side of a shot point, an exponential function is fitted to the amplitude evolution, if at least 4 traces are available. If more than 6 traces exist on the corresponding side, two independent lines are fitted whose results may be interpreted as attenuation near the surface and deeper down. A plot is presented with the fitted logarithm of the amplitudes and the amplitude fit itself. Instead of slope, a Q value is indicated (-1/slope) as well as a r2 value for the ensemble of the two lines.

### **$\textcolor{red}{\text{Picking Menu}}$**

Introductory remarks: try to pick as many of your traces as possible before applying any filter or other data treatment. This is especially true for traces near the shot point, which have in general still high frequencies thus that their waveform is often be strongly affected by frequency or velocity filters.

**$\textcolor{violet}{\text{Manual picks}}$** (keyboard shortcut: **M**): 

Manual picking. For all picks, two clicks with the left mouse button are needed: first you mark the position of the pick, then the uncertainty, which is considered to be symmetric around the pick position, so you may click above or below the pick position. A double click defines as default uncertainty 2 time-samples if no filter was applied to the data and 4 samples if a filter was applied.

To finish, click right mouse button

You may erase picks by clicking the central mouse button (or the mouse wheel). The program searches first the trace nearest to the mouse click and erases then the pick nearest to the click. If no pick is available for the corresponding trace, nothing will happen.

When calling the tool, a stippled line is plotted showing the theoretical air-wave arrivals, which allows avoiding picking this wave as apparent direct wave, when working with near-surface data.

**$\textcolor{violet}{\text{Amplitude picking}}$** (keyboard shortcut: **CTRL-A**): 
Automatic picking tool that is based on finding subsequent maxima for which the relation is maximum after a time that is given by the “maximum velocity” asked for in the dialog box that opens on call. From the maxima found in this way, the program calculates a number of regression lines (number is given by the user) and looks for maximum amplitudes near these regression lines. The half-width of the maximum is a measure of frequency. Signal length allows avoiding the large amplitudes of surface waves from a certain distance on.

**$\textcolor{violet}{\text{Picking STA-LTA}}$** (keyboard shortcut: **SHFT-P**): 
Standard STA-LTA picking. You must give window lengths for STA (short) and LTA (long). You have the possibility to do picking for all traces of the record section or first do a zoom and calculate picks only for the visible traces.

**$\textcolor{violet}{\text{Correlation picking}}$** (keyboard shortcut: **SHFT-C**): 
Do first one manual pick on a good trace which will be used as reference trace. From this pick on, further picks are searched by cross correlation. If picks from other shots exist, they are used by the algorithm as guide (maximum correlation near to other picks made at the same offset). The picks are only calculated for the traces visible on the screen. So, if part of your data is noisy, you may first make a zoom on the good traces, do the correlation picking, then filter the data and do another zoom on the traces not yet picked. Uncertainties are set to 2 samples if no filter was applied to the data, else to 4 samples.

**$\textcolor{violet}{\text{Plot picks}}$** (no keyboard shortcut)
Plots all picks corresponding to actual data gather. Usually, you will not need this tool since picks are plotted automatically.

**$\textcolor{violet}{\text{Move picks}}$** (keyboard shortcut: **CTRL-M**): 
Move picks to other time. Especially useful to adjust automatically done bad picks.

Click left near the pick to be moved. The vertical arrow keyboard keys shift the pick by one sample upward or downward. Arrow key pressed with SHIFT displaces picks by 10 samples, with CTRL, shift picks by 1ms.

When pick is at the correct position, click left on the next pick or use horizontal keyboard keys to go to next trace to the left or to the right.

To finish, click right mouse button

**$\textcolor{violet}{\text{Uncertainty change}}$** (keyboard shortcut: **SHFT-U**): 
Change uncertainties of picks. Use is similar as for moving picks: Vertical keyboard arrows increase or decrease uncertainties. Minimum uncertainty is set to 2 time-samples. Finish with right click on mouse

**$\textcolor{violet}{\text{Erase all picks}}$** (keyboard shortcut: **CTRL-E**): 
Erase all picks of the actually shown seismogram gather (other picks are maintained). Usually used if automatic picking gave too bad results.

**$\textcolor{violet}{\text{Plot calculated times}}$** (no keyboard shortcut)
Activated after using Tomography or if the program at start finds a file called “calc_picks.dat”. The calculated picks are plotted as a connected line. As long as this tool is not pressed again, calculated picks are plotted for all gathers. Deactivation of the tool takes only effect at the next screen refreshment (next gather, next zoom etc.)

**$\textcolor{violet}{\text{Store picks}}$** (keyboard shortcut: **CTRL-S**): 
Store all available picks (not only those of the actual seismogram gather) into two files. The standard file used by PyRefra is “picks.dat”. It has 5 columns:
nr shot point, nr receiver pt, picked time, picked time minus uncertainty, picked time plus uncertainty
nr shot point as in file shots.geo, nr receiver point as in file receivers.geo. The uncertainties are normally considered to be symmetric, but certain automatic picking routines may give different uncertainties before and after the pick.
In the latest version, PyRefra stores automatically all picks each time manual picking is finished with right mouse key. However, automatic picking routines don’t do this yet, mainly because they are too uncertain. It is recommended to use CTRL-S after automatic picking if you are satisfied.

**$\textcolor{violet}{\text{Store Gimli format}}$** (no keyboard shortcut): 
Stores calculated picks in Gimli format (file picks.sgt), used by the Tomography tool. Usually, you do not need this tool, since when calling Tomography, the measured picks are automatically stored in this format before starting inversion.

**ATTENTION**: If your working directory is on an external disk, it seems that losing the contact once is losing it until you restart the program: Therefore, if $\textcolor{violet}{\text{Store Picks}}$ or $\textcolor{violet}{\text{Store GIMLI}}$ detect an error the program is stopped, and you lose only the data not yet stored since the last CTRL-S or manual picking.

### **$\textcolor{red}{\text{Filter Menu}}$**

**$\textcolor{violet}{\text{Filter}}$** (keyboard shortcut: **CTRL-F**): 
Frequency filter. Program calculates average FFT of all shown traces and allows user to define corner frequencies. Click left mouse at position of low-cut frequency, pull line and release at high-cut frequency. If one of the two frequencies is clicked outside the frequency axes (negative frequency or more than maximum frequency), the corresponding fllter is not applied.

After filter frequencies have been defined graphically, a dialog box allows modifying them (e.g., rounding).
If a frequency filter had already been defined before, the dialog box opens immediately, proposing the frequencies of the last applied filter. You may accept, change them manually or ask for graphical redefinition.

**$\textcolor{violet}{\text{Frequency filter for single trace}}$** (no keyboard shortcut): 
Same as standard frequency filter, but you must choose one single trace to filter, by left clicking onto it.

**$\textcolor{violet}{\text{Air wave fk}}$** (keyboard shortcut: **SHFT-A**): 
F-K filtering of all velocities less than 400 m/s.

**$\textcolor{violet}{\text{Velocity filter}}$** (keyboard shortcut: **CTRL-V**): 
F-K filtering choosing interactively the cutting velocity (all smaller velocities will be attenuated).
The cutting velocity is defined interactively by pulling a slider to the desired position. To accept the chosen velocity, release the mouse. The program then asks you whether you want to apply the same filter to all gathers or only to the visible one. The chosen velocity as saved as default velocity for further velocity filtering.

### **$\textcolor{red}{\text{Mute Menu}}$**

**$\textcolor{violet}{\text{Mute/recover trace}}$** (keyboard shortcut: **CTRL-T**): 
Click left on all trace you want to mute (small red star is plotted near the lower end of the trace) or which you want to plot again if it was muted (a small green star is plotted at that place). Finish with right click.

The muting will be applied to all traces having been recorded at the same receiver position.
When the program is finished using “File -> Quit”, a file “receivers modify.dat” is written indicating the muted traces, which may eventually be renamed to “receiver_corrections.dat” to make these changes durable.

**$\textcolor{violet}{\text{Mute air wave}}$** (keyboard shortcut: **A**): 
Mute a stripe around the air wave (mainly thought as preparation before treating data as reflection seismics). In the dialog box give estimated air-wave velocity, the total with of the stripe and the width of the taper to both sides. If you are not happy with the result, go back to original data (SHFT-O) and try again.

**$\textcolor{violet}{\text{Mute before line}}$** (keyboard shortcut: **L**): 
Click left and draw multiple line segments, finish with right click. Data before the traced line are eliminated for each trace before the nearest zero-crossing. Traces that are not crossed by the line (at the right or left side of the screen will not be muted!

**$\textcolor{violet}{\text{Mute after line}}$** (keyboard shortcut: **SHFT-L**): 
Click left and draw multiple line segments, finish with right click. Data after the traced line are eliminated for each trace from the nearest zero-crossing on. Traces that are not crossed by the line (at the right or left side of the screen will not be muted!


