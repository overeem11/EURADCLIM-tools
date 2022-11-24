#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Python script
# Name: VisualizeHDF5ODIMCartopyOSM_GM.py
#
#
## Version 1.0
## Copyright (C) 2022 Aart Overeem
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#
# Description: Script to visualize ODIM HDF5 images (data at a 2-D grid) over Europe.
#              Read ODIM HDF5 files, extract values, read coordinates of center of grid cells. Plot variable on a map using the polygons and associated values.
#              Output is a graphical file with a map. Draws an OpenStreetMap (OSM) or Google Maps (GM) as background. Transparency of polygons can be choosen (alpha), but
#              seems especially useful when zooming in and does not work well at European scale.
#              Note that coastlines and country borders can be plotted on top op the background map. 
#              These are clearly less accurate and can deviate 1-2 km from the ones in the background map.
#              This is not visible on a larger scale, such as Europe.
#              This script is quite generic, since you need to provide the DatasetNr in the ODIM file and you can change label names (variable names), et cetera.
#              So this script can be used to visualize the variable you like from a file which complies to the ODIM-HDF5 standard.
#              Using 300 dpi will result in radar polygons which are a bit distorted ("blurry"), which is especially visible when zooming in.
#              Pdf files lead to very high quality maps, but can be larger than 100 MB, whereas a jpg file at 300 dpi is about 1 MB for example 1 below.
#              Moreover, it takes rather long to produce a pdf output file, whereas producing a jpg output file at 300 dpi takes ~5 seconds.
#              To reduce file size of pdf files after running this script, the program ps2pdf can be used, which converts it to a file with extension ".pdf.pdf".
#              Subsequently, change the extension into ".pdf". Although much reduced, file size will remain relatively large. 
#              A pdf file gives the most accurate results when it comes to plotting radar polygons.
#              Note that the center of pixels (grid cells) are used to plot the polygons. This actually gives a shift of 1 km in latitude and 1 km in longitude, but
#              this is negligible at a map of Europe. Probably it is correct after all, since input coordinates to pcolormesh are interpreted as cell centers.
# Usage: python VisualizeHDF5ODIMCartopyOSM_GM.py [input filename] [output graphical filename, with as extension e.g. "jpg" or "pdf"] [title of map] [legend label text] [color scheme] [scale numbers] [draw country borders] [draw coastlines] [display color for values below lowest number] [colorbar] [extra text] [the resolution of the map in dots per inch] [font size of title] [font size of legend label] [font size of legend] [font size extra text] [draw parallels and meridians] [DatasetNr in ODIM file] [Type of background map: "OSM" for OpenStreetMap & "GM" for Google Maps] [style of Google Maps: e.g. "street" or "satellite"] [value for request: determines resolution of background map / size of text (e.g. city names on map)] [plotting area: minimum longitude, maximum longitude, minimum latitude, maximum latitude] [the alpha blending value, between 0 (transparent) and 1 (opaque) for plotting polygons] [longitude of extra text] [latitude of extra text].
# Example: python VisualizeHDF5ODIMCartopyOSM_GM.py RAD_OPERA_HOURLY_RAINFALL_ACCUMULATION_EURADCLIM/2020/10/RAD_OPERA_HOURLY_RAINFALL_ACCUMULATION_202010021400.h5 figures/RAD_OPERA_HOURLY_RAINFALL_ACCUMULATION_202010021400_EURADCLIM_OSM_France.jpg 'EURADCLIM' '1-h precip depth (mm)' CbF '[5,20,35,50,65,80]' NoDrawCountries NoDrawCoastlines DoNotColorSetUnder ColorBar '(d) 2 Oct 2020 13-14 UTC' 600 29 29 29 26 NoDrawParallelsMeridians '/dataset1' OSM street 12 '[6.5,7.5,43.5,44.2]' 1 6.55 44.15


# Load Python packages:
import sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import h5py
import pyproj
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import pyepsg
import pandas as pd
import cartopy.io.img_tiles as cimgt
import copy



# Parameters from command line:    
InputFileName = sys.argv[1]
OutputFileName = sys.argv[2]
TitlePlot = sys.argv[3]
LabelName = sys.argv[4]                 # Note that this works with mathematical notation as follows (example): 'Radar reflectivity factor (dB$Z$)' (do not use " ").
ScaleType = sys.argv[5]	                # 'CbF' (= Colorblind Friendly) or 'Blues' or 'YellowRed'.
levels = list(map(float, sys.argv[6].strip('[]').split(',')))
LowestValue = float(sys.argv[6].strip('[]').split(',')[0])
DrawCountries = sys.argv[7]		# If DrawCountries is not equal to DrawCountries, the country borders are not drawn.
DrawCoastlines = sys.argv[8]            # If DrawCoastlines is not equal to DrawCoastlines, the coastlines are not drawn.
DoColorSetUnder = sys.argv[9]
ColorBar = sys.argv[10]
ExtraText = sys.argv[11]
dpi = int(sys.argv[12])
FontSizeTitle = sys.argv[13]
FontSizeLegendLabel = sys.argv[14]
FontSizeLegend = sys.argv[15]
FontSizeExtraText = sys.argv[16]
DrawParallelsMeridians = sys.argv[17]   # If DrawParallelsMeridians is not equal to DrawParallelsMeridians, the parallels & meridians are not drawn.
DatasetNr = sys.argv[18]
TypeBackGroundMap = sys.argv[19]
style = sys.argv[20]
ValueRequest = int(sys.argv[21])
extent = list(map(float, sys.argv[22].strip('[]').split(',')))
alpha = float(sys.argv[23])
LonExtraText = float(sys.argv[24])
LatExtraText = float(sys.argv[25])



#############################################################
# 1. Read HDF5 files (ODIM HDF5 format) and obtain polygons.#
#############################################################

DATAFIELD_NAME = DatasetNr + '/data1/data'
FILE_NAME = InputFileName
f = h5py.File(FILE_NAME, mode='r')
# Read metadata:    
Ncols = int(f['/where'].attrs['xsize'])
Nrows = int(f['/where'].attrs['ysize'])

ATTR_NAME = DatasetNr + '/what'
zscale = f[ATTR_NAME].attrs['gain']
zoffset = f[ATTR_NAME].attrs['offset']
nodata = f[ATTR_NAME].attrs['nodata']
undetect = f[ATTR_NAME].attrs['undetect']

# Read data:
dset = f[DATAFIELD_NAME]
RArray = zoffset + zscale * dset[:]

# Read file with coordinates OPERA radar grid:
Grid = np.array(pd.read_csv("CoordinatesHDF5ODIMWGS84.dat", delimiter = " ", dtype="float",header=None))
Xcoor = Grid[:,0]
Ycoor = Grid[:,1]


# Image coordinates to longitude & latitude in degrees:
LonArray = [[np.nan for x in range(Ncols)] for x in range(Nrows)]
LatArray = [[np.nan for x in range(Ncols)] for x in range(Nrows)]
# Obtain image coordinates of surrounding pixels (center of grid cells):
Nrow = Nrows
Ncol = Ncols
for j in range(0,Nrow):
   for i in range(0,Ncol):
      # Coordinate of pixel represents center of pixel.  
      LonArray[j][i] = Xcoor[i+j*Ncols]
      LatArray[j][i] = Ycoor[i+j*Ncols]

# Set data below LowestValue in [scale numbers] to "not available":
RArray[np.isnan(RArray)] = nodata
if DoColorSetUnder!="DoColorSetUnder":
   RArray[RArray < LowestValue] = np.nan
# No data & undetect data values are made "not available".
RArray[RArray == nodata] = np.nan
RArray[RArray == undetect] = np.nan




##################
# 2. Make a plot.#
##################

# Map settings (e.g. projection and extent of area):
plt.rcParams["font.family"] = "serif"
plt.rcParams.update({'font.size': FontSizeLegendLabel}) 
plt.close('all')
transform = ccrs.PlateCarree()
fig = plt.figure(figsize=(16, 8))


# To use OpenStreetMap:
if TypeBackGroundMap=="OSM":
   request = cimgt.OSM()
# To use Google Maps:
if TypeBackGroundMap=="GM":
   request = cimgt.GoogleTiles(style=style)
# Set map:
ax = plt.axes(projection=request.crs)
ax.set_extent(extent)
ax.add_image(request, ValueRequest)

# Plotting it for a chosen projection seems to decrease the graphical quality considerably, so this approach has not been chosen:
#projection = ccrs.epsg(3035)
#ax = plt.axes(projection=projection)



# Choose color scheme for legend:
# 4 classes for 'Blues' does not work well, then you get two white classes! In general, 4 classes seems not to work with get_cmap, so use at least 5 classes.
# Lowest class, e.g. 1 - 26 mm includes the 1 values. Values below 1 are plotted in white or light gray (or another color specified below) depending on the color scheme. 
# The highest value of the highest class is plotted in black or indigo (or another color specified below).
# So for each class its lowest value, i.e. the lowest value at the tick mark, belongs to that class, whereas its highest value belongs to the next class.
# See for colour scales: https://matplotlib.org/2.0.2/examples/color/colormaps_reference.html
if ScaleType=="Blues":
   levels = levels
   cmap = copy.copy(mpl.cm.get_cmap("Blues"))
   colorSetOver = 'indigo'
   colorSetUnder = 'lightgray'
if ScaleType=="YellowRed":
   levels = levels
   cmap = mpl.colors.ListedColormap(['#ffffb2','#fed976','#feb24c','#fd8d3c','#f03b20','#bd0026'])
   colorSetOver = 'black'
   colorSetUnder = 'white'
if ScaleType=="CbF":
   levels = levels
   cmap = mpl.colors.ListedColormap(['#DBEED3','#9CD5C4','#71B5C7','#858AC1','#A2569C','#96344E'])
   colorSetOver = 'black'
   colorSetUnder = 'white'
if ScaleType=="RedBlue":
   levels = levels
   cmap = copy.copy(mpl.cm.get_cmap("seismic_r"))
   colorSetOver = 'black'
   colorSetUnder = 'gray'
if ScaleType=="PiYG":
   levels = levels
   cmap = copy.copy(mpl.cm.get_cmap("PiYG")) 
   colorSetOver = 'darkgreen'   
   colorSetUnder = 'black'
if ScaleType=="PiYG_r":   
   levels = levels
   cmap = copy.copy(mpl.cm.get_cmap("PiYG_r"))
   colorSetOver = 'black'
   colorSetUnder = 'darkgreen'
if ScaleType=="RdYlBu":
   levels = levels
   cmap = copy.copy(mpl.cm.get_cmap("RdYlBu"))
   colorSetOver = 'indigo'
   colorSetUnder = 'black' 
if ScaleType=="Blues_r":
   levels = levels
   cmap = copy.copy(mpl.cm.get_cmap("Blues_r"))
   colorSetOver = 'lightgray'
   colorSetUnder = 'indigo'  


# Plot gridded radar data as colored polygons:
from matplotlib.colors import BoundaryNorm
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
CS3 = plt.pcolormesh(np.asarray(LonArray),np.asarray(LatArray),RArray,cmap=cmap,norm=norm,transform=transform, zorder=2,shading='nearest',alpha=alpha) 
# Set highest class to chosen color "colorSetOver":
CS3.cmap.set_over(colorSetOver,alpha=alpha)



# Plot color bar:
if DoColorSetUnder=="DoColorSetUnder":
   CS3.cmap.set_under(colorSetUnder,alpha=alpha)
   if ColorBar=="ColorBar":
      # Plot color bar:
      font = mpl.font_manager.FontProperties(size=FontSizeLegend)
      cbar = plt.colorbar(pad=0.02,shrink=0.9,extend='both')
      cbar.set_label(LabelName)
      text = ax.yaxis.label
      text.set_font_properties(font)
      cbar.ax.tick_params(labelsize=FontSizeLegend)
else:
   if ColorBar=="ColorBar":
      # Plot color bar:
      font = mpl.font_manager.FontProperties(size=FontSizeLegend)
      cbar = plt.colorbar(pad=0.02,shrink=0.9,extend='max')
      cbar.set_label(LabelName)
      text = ax.yaxis.label
      text.set_font_properties(font)
      cbar.ax.tick_params(labelsize=FontSizeLegend)



# Add natural earth features and borders:
if DrawCountries=="DrawCountries":    
   ax.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=0.3, zorder=2)
if DrawCoastlines=="DrawCoastlines":
    ax.coastlines(resolution='10m', linewidth=0.3, zorder=2)


   
# Plot extra text:
x1, y1 = LonExtraText, LatExtraText
ax.text(x1, y1, ExtraText, color='black', size=FontSizeExtraText, transform=transform)
  


# Draw parallels and meridians:
if DrawParallelsMeridians=="DrawParallelsMeridians":
   gl = ax.gridlines(crs=transform, draw_labels=True, linewidth=1, color='black', linestyle='--')
   gl.top_labels = False
   gl.right_labels = False
   gl.xformatter = LONGITUDE_FORMATTER
   gl.yformatter = LATITUDE_FORMATTER
   gl.xlabel_style = {'size': 15}
   gl.ylabel_style = {'size': 15}



# Plot title:
plt.title(TitlePlot,fontsize=FontSizeTitle)   



# Save figure:
plt.savefig(OutputFileName, bbox_inches = "tight", dpi = dpi)



# Close radar file:
f.close()

