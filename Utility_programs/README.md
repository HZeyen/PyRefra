Python programs useful for preparation and post-treatment of data. For more detailed descriptions see comments within the program files.

**Acquisition-scheme.py** : Plot positions of receivers and shot points and file numbers for acquisition in blocks and shifted geophones (Fig. 6 in publication)

**Extract_velocities.py** : Extract one of three possibilities: the minimum depth of a given velocity along a tomography model; a vertical section; a horizontal section. Program needs library Dialog_HZ.py

**picks_sgt_2_dat.py** : Translates a pick file in gimli *.sgt format to PyRefra format. Input: *.sgt file, receivers.geo and shots.geo. Output: picks.dat.
Coordinates in teh first part of *.sgt file must be the same as those in the *.geo files.

**picks_to_vista-header.py** : Transforms file picks.dat into the format needed for importing picks into VISTA (Schlumberger's reflection seismic data treatment)

**PIK_2_picks_dat.py** : Translates a pick file in Zelt format to PyRefra format

**seg2_renumber.py** : Change number of seg2 files such that numbers are contiguous with the possibility to start at another number than 1

**Tomo_coordinate_projection.py** : Projects local coordinates of a tomography model onto absolute coordinates (e.g., UTM) of the line

**turn_coordinates.py** : Takes files of shot and receiver points in *.geo format with GPS coordinates and creates files receivers.geo and shots.geo 
in a local coordiante system such that the first point is in the local W and the last oint in the local E. In addition, the minimum topography value
is subtracted

**v-tomo2v-RMS.py** :  Transforms a velocity-depth file obtained from tomographic inversion with program PyRefra into RMS velocities vs TWT used with VISTA
(Schlumberger's reflection seismic data treatment)

