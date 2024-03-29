# EURADCLIM-tools: Tools to accumulate and visualize OPERA & EURADCLIM radar data and to perform climatological analyses.

# Introduction
EURADCLIM (EUropean RADar CLIMatology) is a publicly available climatological dataset of 1-h and 24-h precipitation accumulations on a 2-km grid for the period 2013 through 2020. The starting point is the European Meteorological Network (EUMETNET) Operational Program on the Exchange of weather RAdar Information (OPERA) gridded radar dataset of 15-min instantaneous surface rain rates, which is based on data from, on average, 138 ground-based weather radars. After additional removal of non-meteorological echoes by three algorithms, this dataset is merged with precipitation accumulations from potentially ~7700 rain gauges obtained from the European Climate Assessment & Dataset (ECA&D). EURADCLIM covers about 78% of geographical Europe. For an overview of the EURADCLIM project and dataset: https://www.knmi.nl/research/observations-data-technology/projects/euradclim-the-european-climatological-high-resolution-gauge-adjusted-radar-rainfall-dataset.

# Tools
The following tools, written in programming language Python (version 3), are made available in EURADCLIM-tools:
- "AccumulateRadarHDF5ODIMListCount.py": Script to accumulate ODIM HDF5 radar precipitation images where the file names are provided as input or to accumulate ODIM HDF5 radar precipitation images for a given path.
- "ClimatologyRadarHDF5ODIMList.py": Script to perform climatological analysis on ODIM HDF5 radar precipitation images where the file names are provided as input or for a given path. The maximum rainfall, mean rainfall, number of events (= frequency) above a threshold value, or number of events above a threshold value divided by the total number of events (= relative frequency) can be computed.
- "VisualizeHDF5ODIMCartopy.py": Script to visualize ODIM HDF5 images (data at a 2-D grid) over Europe.
- "VisualizeHDF5ODIMCartopyOSM_GM.py": Script to visualize ODIM HDF5 images (data at a 2-D grid) over Europe. Draws an OpenStreetMap or Google Maps as background. Useful when zooming in on a part of Europe.

# Details
Information on usage of these tools is provided in the respective Python script. The provided tools have been tested with the following OPERA composite products:
- Instantaneous Surface Rain Rate
- 1 Hour Rainfall Accumulation

They have also been tested with the 1-h and 24-h precipitation accumulations in EURADCLIM. Note that all these datasets are in the ODIM-HDF5 format on the default OPERA grid of 2 km resolution (Lambert Azimuthal Equal Area projection; 2200 times 1900 grid cells).

Note that some Python libraries need to be installed for these tools to work, but they are not all needed for every tool: cartopy, copy, h5py, matplotlib, natsort, numpy, os, pandas, pyepsg, pyproj, shutil, sys, warnings.

The file "CoordinatesHDF5ODIMWGS84.dat" is needed for the two visualization tools. It contains the coordinates of the center of radar grid cells with longitude (first column) and latitude (second column) in degrees (WGS84). This may also be useful for processing the EURADCLIM radar data for analyses and applications or for other visualization tools.

The file "RAD_OPERA_RAINFALL_RATE_201812110715.h5" is used as a template ODIM-HDF5 file and only needed for tools "AccumulateRadarHDF5ODIMListCount.py" & "ClimatologyRadarHDF5ODIMList.py".

# EURADCLIM datasets
The EURADCLIM datasets of 1-h and 24-h precipitation accumulations are publicly available here: https://doi.org/10.21944/7ypj-wn68 & https://doi.org/10.21944/1a54-gg96. The accompanying scientific article can be found here: https://doi.org/10.5194/essd-15-1441-2023.

# Usage - Example for "VisualizeHDF5ODIMCartopy.py"
```
python VisualizeHDF5ODIMCartopy.py RAD_OPERA_24H_RAINFALL_ACCUMULATION_201305311400.h5 RAD_OPERA_24H_RAINFALL_ACCUMULATION_201305311400_EURADCLIM.jpg 'EURADCLIM' '24-h precip depth (mm)' CbF '[1,10,20,30,40,50,60]' DrawCountries DrawCoastlines DrawLakelines DoNotColorSetUnder ColorBar '(a) End: 31 May 2013 14 UTC' lightgray lightblue 300 24 29 31 21 NA NA NA DrawNorthArrow DrawScaleBar '/dataset1'
```
<img src="RAD_OPERA_24H_RAINFALL_ACCUMULATION_201305311400_EURADCLIM.jpg" alt="drawing" width="500"/>
Figure taken from https://doi.org/10.5194/essd-2022-334. Map made with Natural Earth. Free vector and raster map data &copy naturalearthdata.com.

# Usage - Example for "VisualizeHDF5ODIMCartopyOSM_GM.py"
```
python VisualizeHDF5ODIMCartopyOSM_GM.py RAD_OPERA_HOURLY_RAINFALL_ACCUMULATION_EURADCLIM/2020/10/RAD_OPERA_HOURLY_RAINFALL_ACCUMULATION_202010021400.h5 figures/RAD_OPERA_HOURLY_RAINFALL_ACCUMULATION_202010021400_EURADCLIM_OSM_France.jpg 'EURADCLIM' '1-h precip depth (mm)' CbF '[5,20,35,50,65,80]' NoDrawCountries NoDrawCoastlines DoNotColorSetUnder ColorBar '(d) 2 Oct 2020 13-14 UTC' 600 29 29 29 26 NoDrawParallelsMeridians '/dataset1' OSM street 12 '[6.5,7.5,43.5,44.2]' 1 6.55 44.15
```
<img src="RAD_OPERA_HOURLY_RAINFALL_ACCUMULATION_202010021400_EURADCLIM_OSM_France.jpg" alt="drawing" width="500"/>
Figure taken from https://doi.org/10.5194/essd-2022-334. &copy OpenStreetMap contributors 2022. Distributed under the Open Data Commons Open Database License (ODbL) v1.0.

# Reference
When referring to EURADCLIM-tools, please use:

Aart Overeem. (2022). EURADCLIM-tools (v.1.0). Zenodo. https://doi.org/10.5281/zenodo.7473816
