#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Python script
# Name: AccumulateRadar_netCDF4_CF_Count.py
#
#
## Version 2.0
## Copyright (C) 2026 Aart Overeem
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
# Description: Script to accumulate netCDF4 CF-1.8 radar rainfall images where the file names are provided as input or to accumulate netCDF4 CF-1.8 radar rainfall images for a given path. 
#              Read netCDF4 CF-1.8 radar file, extract rainfall depth, accumulate rainfall depth over provided input files.
#              Output is a netCDF4 CF-1.8 radar file with the accumulated rainfall depth.
#              Use a conversion factor of 1 when input files contain accumulated rainfall (mm).
#              Use a conversion factor of 0.25 when input files contain rainfall intensity in mm per hour.
#              nodata values are set to 0, so all data are accumulated irrespective of their value (e.g. missing data), but
#              the needed minimum availability per radar pixel is employed in the end to decide whether a radar pixel should have data or nodata. 
#              In case of missing files no output file is made. If data criterion is not satisfied, nodata values are present in the output file. 
#              Note that we specifically deal with NaN values by setting them to nodata.

# Usage: python AccumulateRadar_netCDF4_CF.py [output filename] [input file names or path with files] [conversion factor] 
# [needed minimum number of images with data] [path or file names?: choose "path" or "files"] 
# Example (annual precipitation accumulation over the year 2013): python AccumulateRadar_netCDF4_CF.py "RAD_OPERA_HOURLY_RAINFALL_ACCUMULATION_2013_01Y_v4.nc" "/nobackup_2/users/overeem/C3S_EURADCLIM/V4/RAD_OPERA_HOURLY_RAINFALL_ACCUMULATION_EURADCLIM/2013/" 1 7300 path precipitation_amount
# Example (accumulate two files with annual precipitation): python AccumulateRadar_netCDF4_CF.py "RAD_OPERA_HOURLY_RAINFALL_ACCUMULATION_2013_2014_02Y_v4.nc" "RAD_OPERA_HOURLY_RAINFALL_ACCUMULATION_2013_01Y_v4.nc RAD_OPERA_HOURLY_RAINFALL_ACCUMULATION_2014_01Y_v4.nc" 1 2 files precipitation_amount



# Load Python packages:
import sys
import os
import numpy as np
import shutil
from pathlib import Path
import warnings
#warnings.filterwarnings("ignore")
import natsort
import xarray as xr
import netCDF4 as nc


# Parameters from command line:    
OutputFileName = sys.argv[1]
InputFileNames = sys.argv[2]
ConversionFactor = float(sys.argv[3])
MinImages = int(sys.argv[4])
PathOrFileNames = sys.argv[5]
DatasetVariableName = sys.argv[6]


if PathOrFileNames=="files":
    # Split string with input file names in string with multiple lines:
    pathlist_temp = InputFileNames.split()
    # Sort file names:
    pathlist = natsort.natsorted(pathlist_temp)

if PathOrFileNames=="path":
    pathlist_temp = Path(InputFileNames).glob('**/*.nc')
    # Sort file names:
    pathlist = natsort.natsorted(pathlist_temp)

if PathOrFileNames!="files" and PathOrFileNames!="path":
    print("Accumulation cannot be performed. Please specify whether a list of files (files) or a directory path with files (path) is supplied!")
    sys.exit(0)



########################################################################
#  Read netCDF4 CF-1.8 radar files and accumulate precipitation amount.#
########################################################################


i = 0
radardata = radardata_temp = Count = []
for path in pathlist:
    if os.path.getsize(path) > 0:               # Is used to check the size of specified path. It returns the size of specified path in bytes.
    						# The method raises OSError if the file does not exist or is somehow inaccessible.
       i = i + 1
       print(path)
       if i==1:
          # Open file:
          f = xr.open_dataset(path)
          InputFileName = str(path)
          # Read data:
          nodata = f[DatasetVariableName].encoding.get("_FillValue")
          radardata = f[DatasetVariableName]
          radardata = radardata.squeeze(drop=True) 
          # Replace nan with nodata, and set nodata to 0:
          radardata = radardata.fillna(nodata)
          radardata = xr.where(radardata == nodata, 0, radardata)
          # Count for each radar pixel the number of valid values, i.e. not equal to nodata or nan:
          Count = f[DatasetVariableName]
          Count = Count.squeeze(drop=True) 
          Count = Count.fillna(nodata)
          Count = xr.where(Count == nodata, 0, 1)
          f.close()
       else:
          # Open file:
          f = xr.open_dataset(path)
          InputFileName = str(path)
          # Read data:  
          nodata = f[DatasetVariableName].encoding.get("_FillValue")
          radardata_temp = []  
          radardata_temp = f[DatasetVariableName]
          radardata_temp = radardata_temp.squeeze(drop=True) 
          # Replace nan with nodata, and set nodata to 0:
          radardata_temp = radardata_temp.fillna(nodata)
          radardata_temp = xr.where(radardata_temp==nodata, 0, radardata_temp)
          radardata = xr.ufuncs.add(radardata, radardata_temp)
          # Count for each radar pixel the number of valid values, i.e. not equal to nodata or nan:
          Count_temp = f[DatasetVariableName]
          Count_temp = Count_temp.squeeze(drop=True) 
          Count_temp = Count_temp.fillna(nodata)
          Count = Count + xr.where(Count_temp!=nodata, 1, 0)
          f.close()



###################################################
# Make output radar file in netCDF4 CF-1.8 format.#
###################################################
if i > 0:        # Implies that at least 1 radar image can be read, but can contain nodata only.
   # Apply conversion factor:
   radardata = radardata * ConversionFactor
   # Apply data availability criterion separately per pixel. If availability is too low, set radardata to nodata:
   radardata = xr.where(Count>=MinImages, radardata, nodata)
   # Remove output file:
   try:
       os.remove(OutputFileName)
   except OSError:
       pass
   # Copy input file to netCDF4 output file:       
   shutil.copy(InputFileName, OutputFileName)  
   # Note that metadata is copied from the last file.   
   
   if hasattr(radardata, 'values'):
       radardata_array = radardata.values
   else:
       radardata_array = np.array(radardata)

   # Open netCDF file in writing mode and replace the values for DatasetVariableName with the computed values:
   with nc.Dataset(OutputFileName, 'r+') as nc_file:
       var = nc_file.variables[DatasetVariableName]
       var[:] = radardata_array

   # Because the netCDF4 CF-1.8 file only contains "time", no distinction between start or end date & time can be made. Hence, "time" is just the time from the file that has been copied, i.e.,
   # the last file that could be read.
   # In addition, note that the netCDF4 CF-1.8 files only provide the end date & time of observation. So, "time" actually represents the end date & time of the last file which could be read.


if i==0:
   print("All radar files are empty. Hence, an output radar file could not be produced.")



