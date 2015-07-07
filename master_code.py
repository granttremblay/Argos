# Master code for the Summer project work of 2015
# Launch CIAO before commencing

# Imports

import numpy
import os
import sys
import glob
import pyfits

# Directory specification

directory = raw_input("Enter directory for work: ")
os.system("mkdir "+directory)
os.system("cd "+directory)

# Did this part to try to avoid invalid entries

while True:
	try:
		obsID_choice = raw_input("Enter obsID's: ")
		os.system("download_chandra_obsid "+obsID_choice)
	except ValueError:
		print "Enter a valid obsID"
	else:
		break

# chandra_repro and flux_obs steps

def lvl2_data():
	os.system("chandra_repro indir=%s outdir=%s/repro verbose=5 cleanup=no" % )

# Aligning observations to WCS
os.system("mkdir reprojected_data")
os.system("reproject_obs '*/repro/*evt2.fits' reprojected_data/ verbose=5")

# Create exposure corrected mosaic
os.system("mkdir exposure_corrected_mosaic")
os.system("flux_obs 'reprojected_data/*reproj_evt.fits' exposure_corrected_mosaic/ bands=broad,csc bin=1 verbose=5")

# Create the obsID filenames
os.system("cd reprojected_data/")
os.system("vim main_filenames.txt")
