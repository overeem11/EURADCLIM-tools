#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Python script
# Name: VisualizeHDF5ODIMCartopy.py
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
#              Output is a graphical file with a map.
#              This script is quite generic, since you need to provide the DatasetNr in the ODIM file and you can change label names (variable names), et cetera.
#              So this script can be used to visualize the variable you like from a file which complies to the ODIM-HDF5 standard.
#              Using 300 dpi will result in radar polygons which are a bit distorted ("blurry"), which is especially visible when zooming in.
#              Pdf files lead to very high quality maps, but can be larger than 100 MB, whereas a jpg file at 300 dpi is about 1 MB for example 1 below.
#              Moreover, it takes rather long to produce a pdf output file, whereas producing a jpg output file at 300 dpi takes ~10 seconds.
#              To reduce file size of pdf files after running this script, the program ps2pdf can be used, which converts it to a file with extension ".pdf.pdf".
#              Subsequently, change the extension into ".pdf". Although much reduced, file size will remain relatively large. 
#              A pdf file gives the most accurate results when it comes to plotting radar polygons.
#              Note that the center of pixels (grid cells) are used to plot the polygons. This actually gives a shift of 1 km in latitude and 1 km in longitude, but
#              this is negligible at a map of Europe. Probably it is correct after all, since input coordinates to pcolormesh are interpreted as cell centers.
# Usage: python VisualizeHDF5ODIMCartopy.py [input filename] [output graphical filename, with as extension e.g. "jpg" or "pdf"] [title of map] [legend label text] [color scheme] [scale numbers] [draw country borders] [draw coastlines] [draw lake lines] [display color for values below lowest number] [colorbar] [extra text] [color of land] [color of ocean, rivers & lakes] [the resolution of the map in dots per inch] [font size of title] [font size of legend label] [font size of legend] [font size extra text] [draw parallels and meridians] [draw rivers] [draw provinces] [draw north arrow] [plot scale bar 500 km] [DatasetNr in ODIM file]
# Example: python VisualizeHDF5ODIMCartopy.py RAD_OPERA_24H_RAINFALL_ACCUMULATION_EURADCLIM/2013/05/RAD_OPERA_24H_RAINFALL_ACCUMULATION_201305311400.h5 RAD_OPERA_24H_RAINFALL_ACCUMULATION_201305311400_EURADCLIM.jpg 'EURADCLIM' '24-h precip depth (mm)' CbF '[1,10,20,30,40,50,60]' DrawCountries DrawCoastlines DrawLakelines DoNotColorSetUnder ColorBar '(a) End: 31 May 2013 14 UTC' lightgray lightblue 300 24 29 31 21 NA NA NA DrawNorthArrow DrawScaleBar '/dataset1'


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
DrawLakelines = sys.argv[9]             # If DrawLakelines is not equal to DrawLakelines, the lake lines are not drawn.
DoColorSetUnder = sys.argv[10]
ColorBar = sys.argv[11]
ExtraText = sys.argv[12]
ColorLand = sys.argv[13]
ColorOceanRiverLakes = sys.argv[14]
dpi = int(sys.argv[15])
FontSizeTitle = sys.argv[16]
FontSizeLegendLabel = sys.argv[17]
FontSizeLegend = sys.argv[18]
FontSizeExtraText = sys.argv[19]
DrawParallelsMeridians = sys.argv[20]   # If DrawParallelsMeridians is not equal to DrawParallelsMeridians, the parallels & meridians are not drawn.
DrawRivers = sys.argv[21]               # If DrawRivers is not equal to DrawRivers, the rivers are not drawn.
DrawProvinces = sys.argv[22]            # If DrawProvinces is not equal to DrawProvinces, the departments/provinces are not drawn.
DrawNorthArrow = sys.argv[23]           # If DrawNorthArrow is not equal to NorthArrow, the north arrow is not drawn.
DrawScaleBar = sys.argv[24]             # If DrawScaleBar is not equal to DrawScaleBar, the scale bar are not drawn.
DatasetNr = sys.argv[25]
ColorExtraText = 'black'
if len(sys.argv)==27:
   ColorExtraText = sys.argv[26]



# Function for plotting scale bar (taken from https://stackoverflow.com/questions/32333870/how-can-i-show-a-km-ruler-on-a-cartopy-matplotlib-plot)
def scale_bar(ax, length=None, location=(0.44, 0.91), linewidth=3):
    """
    ax is the axes to draw the scalebar on.
    length is the length of the scalebar in km.
    location is center of the scalebar in axis coordinates.
    (ie. 0.5 is the middle of the plot)
    linewidth is the thickness of the scalebar.
    """
    #Get the limits of the axis in lat long
    llx0, llx1, lly0, lly1 = ax.get_extent(ccrs.PlateCarree())
    #Make tmc horizontally centred on the middle of the map,
    #vertically at scale bar location
    sbllx = (llx1 + llx0) / 2
    sblly = lly0 + (lly1 - lly0) * location[1]
    tmc = ccrs.TransverseMercator(sbllx, sblly, approx=True)  # Aart Overeem (KNMI): added ", approx=True" because of warning.
    #Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(tmc)
    #Turn the specified scalebar location into coordinates in metres
    sbx = x0 + (x1 - x0) * location[0]
    sby = y0 + (y1 - y0) * location[1]

    #Calculate a scale bar length if none has been given
    #(Theres probably a more pythonic way of rounding the number but this works)
    if not length: 
        length = (x1 - x0) / 5000 #in km
        ndim = int(np.floor(np.log10(length))) #number of digits in number
        length = round(length, -ndim) #round to 1sf
        #Returns numbers starting with the list
        def scale_number(x):
            if str(x)[0] in ['1', '2', '5']: return int(x)        
            else: return scale_number(x - 10 ** ndim)
        length = scale_number(length) 

    #Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbx - length * 500, sbx + length * 500]
    #Plot the scalebar
    ax.plot(bar_xs, [sby, sby], transform=tmc, color='k', linewidth=linewidth)
    #Plot the scalebar label
    ax.text(sbx, sby, str(length) + ' km', transform=tmc,
            horizontalalignment='center', verticalalignment='bottom', size=19)



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
projection = ccrs.epsg(3035)     #epsg:3035 = ETRS89 / ETRS-LAEA
fig = plt.figure(figsize=(8, 8))
extent = [-10, 30, 32.7, 73]
ax = plt.axes(projection=projection)
ax.set_extent(extent)



# Choose color scheme for legend:
# 4 classes for 'Blues' does not work well, then you get two white classes! In general, 4 classes seems not to work with get_cmap, so use at least 5 classes.
# Lowest class, e.g. 1 - 26 mm includes the 1 values. Values below 1 are plotted in white or light gray depending on the color scheme. 
# The highest value of the highest class is plotted in black or indigo.
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
if ScaleType=="Blues_r":
   levels = levels
   cmap = copy.copy(mpl.cm.get_cmap("Blues_r"))
   colorSetOver = 'tomato'
   colorSetUnder = 'indigo'     



# Plot gridded radar data as colored polygons:
from matplotlib.colors import BoundaryNorm
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
CS3 = plt.pcolormesh(np.asarray(LonArray),np.asarray(LatArray),RArray,cmap=cmap,norm=norm,transform=transform, zorder=2,shading='nearest') 
# Set highest class to chosen color "colorSetOver":
CS3.cmap.set_over(colorSetOver)



# Plot color bar:
if DoColorSetUnder=="DoColorSetUnder":
   CS3.cmap.set_under(colorSetUnder)
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
ax.add_feature(cartopy.feature.LAND, facecolor=ColorLand)
ax.add_feature(cartopy.feature.OCEAN, facecolor=ColorOceanRiverLakes)
ax.add_feature(cartopy.feature.LAKES, facecolor=ColorOceanRiverLakes, linewidth=0.00001,zorder=1)
if DrawRivers=="DrawRivers":
   ax.add_feature(cartopy.feature.RIVERS, edgecolor=ColorOceanRiverLakes, linewidth=1.2, zorder=2)
if DrawProvinces=="DrawProvinces":
   ax.add_feature(cartopy.feature.STATES.with_scale('10m'), linewidth=0.3, zorder=2, edgecolor='gray')
if DrawCountries=="DrawCountries":    
   ax.add_feature(cartopy.feature.BORDERS, linestyle='-', linewidth=0.3, zorder=2)
if DrawCoastlines=="DrawCoastlines":
    ax.coastlines(resolution='10m', linewidth=0.3, zorder=2)
if DrawLakelines=="DrawLakelines":
   ax.add_feature(cartopy.feature.LAKES, edgecolor='black', linewidth=0.3, facecolor="none",zorder=2)


   
# Plot extra text:
x1, y1 = -9, 33.5
ax.text(x1, y1, ExtraText, color=ColorExtraText, size=FontSizeExtraText, transform=transform)
  


# Draw parallels and meridians:
if DrawParallelsMeridians=="DrawParallelsMeridians":
   gl = ax.gridlines(crs=transform, draw_labels=True, linewidth=1, color='black', linestyle='--')
   gl.top_labels = False
   gl.right_labels = False
   gl.xformatter = LONGITUDE_FORMATTER
   gl.yformatter = LATITUDE_FORMATTER
   gl.xlabel_style = {'size': 15}
   gl.ylabel_style = {'size': 15}



# Plot North arrow: 
if DrawNorthArrow=="DrawNorthArrow":
   ax.annotate('', xy=(27.9, 37.3), xytext=(27.9, 33.3),xycoords=transform._as_mpl_transform(ax),size=25,arrowprops=dict(facecolor="black",arrowstyle="fancy"))
   ax.text(27.95, 37.3, 'N', fontsize=31, transform=transform)



# Plot title:
plt.title(TitlePlot,fontsize=FontSizeTitle)   



# Plot scale bar:
if DrawScaleBar=="DrawScaleBar":
   scale_bar(ax, 500)



# Save figure:
plt.savefig(OutputFileName, bbox_inches = "tight", dpi = dpi)



# Close radar file:
f.close()

