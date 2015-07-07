#!/usr/bin/env python
# Master code for the Summer project work of 2015
# Remember to launch CIAO before commencing
# I am replacing most of the CIAO code with print statements for now, later find and replace print with os.system.  Also only working for one obsID.

# Imports

import numpy
import os
import sys
import glob
import pyfits
import subprocess
import os.path

# Directory specification

directory = raw_input("Enter directory for work: ")
os.system("mkdir "+directory)
os.chdir(directory)

# Loop for finding obsID's (These ValueError's need to be changed to CIAO errors #checker)

while True:
	try:
		os.system("find_chandra_obsID "+raw_input("Enter object name (make sure it is in quotes): "))
		print "Found the object obsID's"
	except ValueError:
		print "Enter a valid object ID"
	else:
		break

# Loop for obsID entry and Chandra download

while True:
	try:
		obsID_choice = raw_input("Enter obsID's: ")
		print "download_chandra_obsid "+obsID_choice
	except ValueError:
		print "Enter a valid obsID"
	else:
		break

# Reduce data to lvl 2 #checker

print ("chandra_repro indir='{0}' outdir='{0}'/repro verbose=5 cleanup=no" .format(obsID_choice))

# Align observations to the same WCS space

os.system("mkdir reprojected_data")
print "reproject_obs '*/repro/*evt2.fits' reprojected_data/ verbose=5"

# Create exposure corrected merged mosaic

os.system("mkdir exposure_corrected_mosaic")
print "flux_obs 'reprojected_data/*reproj_evt.fits' exposure_corrected_mosaic/ bands=broad,csc bin=1 verbose=5"

# Opening broad_flux.img

os.chdir("exposure_corrected_mosaic")
print "ds9 broad_flux.img"  # Opens broad_flux.img in DS9 for making regions
print "Create circular regions around the target and around a background portion of the image in DS9.  Save the target region as cluster.reg and the background region as background.reg.  Save both in /reprojected_data.  Make sure to save the regions in ciao format."

# Waiting for the files to be saved.  The user must indicate when to continue.

raw_input("Press Enter to continue...")
os.chdir("../reprojected_data")

# Checking for cluster.reg and background.reg (#checker)

# while True:
#	try:
#		# os.path.isfile("cluster.reg")
#	except:
#		print "The file cluster.reg was not found in reprojected_data/"
#	else:
#		break

# Create the .txt files with the filenames

line1 =
line2 =

obj_file = open("object_filenames.txt", 'a')
obj_file.write("line1 \n line2 \n line3")
obj_file.close()

