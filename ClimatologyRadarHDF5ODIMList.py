#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Python script
# Name: ClimatologyRadarHDF5ODIMList.py
#
#
## Version 1.1
## Copyright (C) 2025 Aart Overeem
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
# Description: Script to perform climatological analysis on ODIM HDF5 radar rainfall images where the file names are provided as input or for a given path. 
#              Read ODIM HDF5 radar file, extract rainfall depth, perform analysis over provided input files.
#              Output is an ODIM HDF5 radar file with the output variable. The field "dataset2", the quality indicator QIND or count field, is removed.
#              Use a conversion factor of 1 when input files contain accumulated rainfall (mm).
#              Use a conversion factor of 0.25 when input files contain rainfall intensity in mm per hour.
#              Note that conversion factor is only needed for "max" and "mean". In case of computation of (relative) frequencies,
#              choose the threshold value in the unit of the original data, e.g. millimeters in case of accumulations and millimeters per hours in case of intensities.
#              nan, nodata and undetect are set to 0, so all data are used irrespective of their value (e.g. missing data), but
#              the needed minimum availability per radar pixel is employed in the end to decide whether a radar pixel should have data or nodata.  
#              The following output variable, i.e. analysis, can be chosen to be performed over the list of files / files in a given path: 
#              maximum rainfall ("max"), mean rainfall ("mean"), number of events (= frequency) above a threshold value ("freq"), number of events above a threshold value
#              divided by the total number of events (= relative frequency; "relfreq").
#              Also in case of missing files an output file is made. If data criterion is not satisfied nodata values are present in the output file. 
#              Note that we specifically deal with NaN values by setting them to nodata.
#              Note that this script assumes the default gain = 1.0 & offset = 0.0.
# Usage: python ClimatologyRadarHDF5ODIMList.py [output filename] [input file names or path with files] [end date for metadata what] [end time for metadata what] [conversion factor] 
# [needed minimum number of images with data] [input file which is used to construct output file in case all input files are not valid] [path or file names?: choose "path" or "files"] 
# [operation: max/mean/freq/relfreq] [threshold value (only needs to be provided when operation is "freq" or "relfreq")]
# Example (mean hourly rainfall over 2013-2022): python ClimatologyRadarHDF5ODIMList.py "RAD_OPERA_HOURLY_RAINFALL_ACCUMULATION_2013_2022_01Y_Mean_EURADCLIM.h5" "RAD_OPERA_HOURLY_RAINFALL_ACCUMULATION_EURADCLIM" "20221231" "230000" 1 73040 RAD_OPERA_RAINFALL_RATE_201812110715.h5 path mean
# Example (relative frequency of exceedance of 5 mm in an hour over 2013-2022): python ClimatologyRadarHDF5ODIMList.py "RAD_OPERA_HOURLY_RAINFALL_ACCUMULATION_2013_2022_01Y_RelFreq_5mm_EURADCLIM.h5" "RAD_OPERA_HOURLY_RAINFALL_ACCUMULATION_EURADCLIM" "20221231" "230000" 1 73040 RAD_OPERA_RAINFALL_RATE_201812110715.h5 path relfreq 5
# Example (frequency of exceedance of 25 mm in an hour over 2013-2022): python ClimatologyRadarHDF5ODIMList.py "RAD_OPERA_HOURLY_RAINFALL_ACCUMULATION_2013_2022_01Y_Freq_EURADCLIM.h5" "RAD_OPERA_HOURLY_RAINFALL_ACCUMULATION_EURADCLIM" "20221231" "230000" 1 73040 RAD_OPERA_RAINFALL_RATE_201812110715.h5 path freq 25 
# Example (maximum 24-h rainfall over 2013-2022): python ClimatologyRadarHDF5ODIMList.py "RAD_OPERA_24H_RAINFALL_ACCUMULATION_2013_2022_01Y_Max_EURADCLIM.h5" "RAD_OPERA_24H_RAINFALL_ACCUMULATION_EURADCLIM" "20221231" "230000" 1 73040 RAD_OPERA_RAINFALL_RATE_201812110715.h5 path max


# Load Python packages:
import sys
import os
import numpy as np
np.seterr(divide='ignore', invalid='ignore')    #  To handle "divide by zero" or "divide by NaN" and ignore it.
import shutil
from pathlib import Path
import warnings
#warnings.filterwarnings("ignore")
import h5py
import natsort


# Parameters from command line:    
OutputFileName = sys.argv[1]
InputFileNames = sys.argv[2]
Date = sys.argv[3]
Time = sys.argv[4]
ConversionFactor = float(sys.argv[5])
MinImages = int(sys.argv[6])
InputFileNameNodata = sys.argv[7]
PathOrFileNames = sys.argv[8]
Operation = sys.argv[9]
if Operation=="freq" or Operation=="relfreq":
   threshold = float(sys.argv[10])

if Operation!="freq" and Operation!="relfreq" and Operation!="max" and Operation!="mean":
    print("Program cannot be used. Please specify which operation should be performed (max/mean/freq/relfreq)!")
    sys.exit(0)

if PathOrFileNames=="files":
    # Split string with input file names in string with multiple lines:
    pathlist_temp = InputFileNames.split()
    # Sort file names:
    pathlist = natsort.natsorted(pathlist_temp)

if PathOrFileNames=="path":
    pathlist_temp = Path(InputFileNames).glob('**/*.h5')
    # Sort file names:
    pathlist = natsort.natsorted(pathlist_temp)

if PathOrFileNames!="files" and PathOrFileNames!="path":
    print("Program cannot be used. Please specify whether a list of files (files) or a directory path with files (path) is supplied!")
    sys.exit(0)



################################################################################
# Read HDF5 radar files (ODIM HDF5 format) and perform climatological analysis.#
################################################################################


i = 0
DATAFIELD_NAME = '/dataset1/data1/data'
radardata = radardata_temp = Count = []
for path in pathlist:
    if os.path.getsize(path) > 0:               # Is used to check the size of specified path. It returns the size of specified path in bytes.
    						# The method raises OSError if the file does not exist or is somehow inaccessible.
       if h5py.is_hdf5(path):                   # Check that a file is a valid HDF5 file.
          i = i + 1
          print(path)
          if i==1:
             # Open file:
             f = h5py.File(path, "r")
             InputFileName = str(path)
             # Read data:
             nodata = f['/dataset1/what'].attrs['nodata']
             undetect = f['/dataset1/what'].attrs['undetect']
             radardata = f[DATAFIELD_NAME][()]
             truth_table = radardata==undetect
             indices = np.where(truth_table)
             radardata[indices] = 0
             # Replace nan with nodata:
             radardata[np.isnan(radardata)] = nodata
             truth_table = radardata==nodata
             indices = np.where(truth_table)
             radardata[indices] = 0
             startdate = f['/dataset1/what'].attrs['startdate']
             starttime = f['/dataset1/what'].attrs['starttime']
             #np.savetxt("array.txt",radardata, fmt="%s")
             # Count for each radar pixel the number of valid values, i.e. not equal to nodata or nan:
             Count = f[DATAFIELD_NAME][()]
             Count[np.isnan(Count)] = nodata
             truth_table = Count!=nodata
             indices = np.where(truth_table)
             Count[indices] = 1
             truth_table = Count==nodata
             indices = np.where(truth_table)
             Count[indices] = 0
             if Operation=="freq" or Operation=="relfreq":
                Count_freq = f[DATAFIELD_NAME][()]
                Count_freq[np.isnan(Count_freq)] = nodata
                truth_table = Count_freq!=nodata
                indices = np.where(truth_table)
                Count_freq[indices] = 0
                truth_table = Count_freq==nodata
                indices = np.where(truth_table)
                Count_freq[indices] = 0
                truth_table_freq = radardata > threshold
                indices_freq = np.where(truth_table_freq)
                Count_freq[indices_freq] = 1 
             f.close()
          else:
             # Open file:
             f = h5py.File(path, mode='r')
             # Read data:  
             nodata = f['/dataset1/what'].attrs['nodata']
             undetect = f['/dataset1/what'].attrs['undetect']
             radardata_temp = []  
             radardata_temp = f[DATAFIELD_NAME][()]
             truth_table = radardata_temp==undetect
             indices = np.where(truth_table)
             radardata_temp[indices] = 0
             # Replace nan with nodata:
             radardata_temp[np.isnan(radardata_temp)] = nodata
             truth_table = radardata_temp==nodata
             indices = np.where(truth_table)
             radardata_temp[indices] = 0
             if Operation=="max":
                radardata = np.maximum(radardata,radardata_temp)
             if Operation=="mean":
                radardata = np.add(radardata,radardata_temp)
             if Operation=="freq" or Operation=="relfreq":
                truth_table_freq = radardata_temp > threshold
                indices_freq = np.where(truth_table_freq)
                Count_freq[indices_freq] = Count_freq[indices_freq] + 1   
             enddate = f['/dataset1/what'].attrs['enddate']
             endtime = f['/dataset1/what'].attrs['endtime']
             # Count for each radar pixel the number of valid values, i.e. not equal to nodata or nan:
             Count_temp = f[DATAFIELD_NAME][()]
             Count_temp[np.isnan(Count_temp)] = nodata
             truth_table = Count_temp!=nodata
             indices = np.where(truth_table)
             Count[indices] = Count[indices] + 1
             f.close()



##############################################
# Make output radar file in ODIM HDF5 format.#
##############################################
if i > 0:        # Implies that at least 1 radar image can be read, but can contain nodata only.
   if Operation=="mean" or Operation=="max":
      # Apply conversion factor:
      radardata = radardata * ConversionFactor
   # Compute mean rainfall depth:
   if Operation=="mean":
      radardata = radardata / Count
   # Compute relative frequency of exceeding a threshold value:
   if Operation=="relfreq":
      radardata = []
      radardata = Count_freq / Count
   if Operation=="freq":
      radardata = []
      radardata = Count_freq
   # Apply data availability criterion separately per pixel. If availability is too low, set radardata to nodata:
   truth_table = Count<MinImages
   indices = np.where(truth_table)    
   radardata[indices] = nodata
   # Remove output file:
   try:
       os.remove(OutputFileName)
   except OSError:
       pass
   # Copy input file to HDF5 output file:       
   shutil.copy(InputFileName, OutputFileName)  
   # Note that metadata is copied from the first file and may hence not be entirely representative for the entire accumulation period (e.g. "how" (which radars)).   
   # Make output file:
   hf = h5py.File(OutputFileName, "a")
   # Remove data field:	
   del hf[DATAFIELD_NAME]
   # Remove quality indicator / count field if it exists:
   DATAFIELD_NAME_2 = '/dataset2'
   NodeExists = DATAFIELD_NAME_2 in hf
   if NodeExists==True:
     del hf[DATAFIELD_NAME_2]
   # Create data field including attributes and write to HDF5 output file:
   dset = hf.create_dataset(DATAFIELD_NAME, data=radardata, compression="gzip", compression_opts=6)
   dset.attrs["CLASS"] = np.string_("IMAGE ")
   dset.attrs["IMAGE_VERSION"] = np.string_("1.2 ") 
   # Modify attributes:
   dset = hf["dataset1/what"]
   # If all files are available the start and end date & time will be correct.
   # Note that for the start date and time the first file which could be read is used, whereas for the end date and time the last file which could be read is used.
   # So in case of missing files, the start and/or end date & time can become different. In case of a missing file between start and end date & time, this cannot be seen
   # in the attributes enddate, endtime, startdate, and starttime.
   if i > 1:
       dset.attrs.modify('enddate',enddate)
       dset.attrs.modify('endtime',endtime)
   if i==1:
       dset.attrs.modify('enddate',np.string_(Date))
       dset.attrs.modify('endtime',np.string_(Time))
   dset.attrs.modify('startdate',startdate)
   dset.attrs.modify('starttime',starttime)
   dset.attrs.modify('quantity',np.string_("CLIM "))
   #
   # Change date and time in general attributes:
   dset = hf["what"]
   dset.attrs.modify('date',np.string_(Date))
   dset.attrs.modify('time',np.string_(Time))
   hf.close()


if i==0:                 # If none of the input files can be read, make an output file with metadata and the radardata field being entirely nodata:
   print("All radar files are empty. Hence, an output radar image could not be produced, but a file with nodata entries is produced.")
   if os.path.getsize(InputFileNameNodata) > 0:               # Is used to check the size of specified path. It returns the size of specified path in bytes.
                                                              # The method raises OSError if the file does not exist or is somehow inaccessible.
      if h5py.is_hdf5(InputFileNameNodata):                   # Check that a file is a valid HDF5 file.
         # Open file:
         f = h5py.File(InputFileNameNodata, "r")
         # Read data:
         nodata = f['/dataset1/what'].attrs['nodata']
         radardata = f[DATAFIELD_NAME][()]
         truth_table = radardata!=nodata
         indices = np.where(truth_table)
         radardata[indices] = nodata
         f.close()
         # Remove output file:
         try:
             os.remove(OutputFileName)
         except OSError:
             pass
         # Copy input file to HDF5 output file:       
         shutil.copy(InputFileNameNodata, OutputFileName) 
         # Note that metadata is copied from a file from another period and may hence not be (entirely) representative for the entire accumulation period (e.g. "how" (which radars)).   
         # Make output file:
         hf = h5py.File(OutputFileName, "a")
         # Remove data field:	
         del hf[DATAFIELD_NAME]
         # Remove quality indicator / count field if it exists:
         DATAFIELD_NAME_2 = '/dataset2'
         NodeExists = DATAFIELD_NAME_2 in hf
         if NodeExists==True:
           del hf[DATAFIELD_NAME_2]
         # Create data field including attributes and write to HDF5 output file:
         dset = hf.create_dataset(DATAFIELD_NAME, data=radardata, compression="gzip", compression_opts=6)
         dset.attrs["CLASS"] = np.string_("IMAGE ")
         dset.attrs["IMAGE_VERSION"] = np.string_("1.2 ") 
         # Modify attributes:
         dset = hf["dataset1/what"]
         # It seems that the Time is usually e.g. 220500 instead of 220000. Here, we just use 220000.
         # Moreover, we fill in the end time of observation for both the end and start time.
         # Since these metadata are generally not used and the end data and time in the general node "what" is OK, we decide to leave it like this.
         # Note that this only occurs if all files for a given accumulation interval cannot be opened.
         dset.attrs.modify('enddate',np.string_(Date))
         dset.attrs.modify('endtime',np.string_(Time))
         dset.attrs.modify('startdate',np.string_(Date))
         dset.attrs.modify('starttime',np.string_(Time))
         dset.attrs.modify('quantity',np.string_("CLIM "))
         #
         # Change date and time in general attributes:         
         dset = hf["what"]
         # This is correct, because this is the end date and time:
         dset.attrs.modify('date',np.string_(Date))
         dset.attrs.modify('time',np.string_(Time))
         hf.close()
      else:
         print("Nodata input file cannot be opened! No output file is constructed!")
   else:
      print("Nodata input file cannot be opened! No output file is constructed!")






