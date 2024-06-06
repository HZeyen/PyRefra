Python programs useful for preparation and post-treatment of data

**Extract_velocity_isoline.py** : Extract the minimum depth of a given velocity along a tomography model

**picks_sgt_2_dat.py** : Translates a pick file in gimli *.sgt format to PyRefra format. Input: *.sgt file, receivers.geo and shots.geo. Output: picks.dat.
Coordinates in teh first part of *.sgt file must be the same as those in the *.geo files.

**picks_to_vista-header.py** : Transforms file picks.dat into the format needed for importing picks into VISTA (Schlumberger's reflection seismic data treatment)

**PIK_2_picks_dat.py** : Translates a pick file in Zelt format to PyRefra format

**seg2_renumber.py** : Change number of seg2 files such that numbers are contiguous with th possibility to start at another number than 1

**Tomo_coordinate_projection.py** : Projects local coordinates of a tomography model onto absolute coordinates (e.g., UTM) of the line

**turn_coordinates.py** : Takes files of shot and receiver points in *.geo format with GPS coordinates and creates files receivers.geo and shots.geo 
in a local coordiante system such that the first point is in the local W and the last oint in the local E. In addition, the minimum topography value
is subtracted

**v-tomo2v-RMS.py** :  Transforms a velocity-depth file obtained from tomographic inversion with program PyRefra into RMS velocities vs TWT used with VISTA
(Schlumberger's reflection seismic data treatment)

