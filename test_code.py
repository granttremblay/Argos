#!/usr/bin/env python
# File for testing code segments

import numpy
import os
import sys
import glob
import pyfits
import subprocess
from astropy.io import fits

obsID_li = []

def obsID_selection():

	count = 1
	while 1:
		obsID = raw_input("Enter obsID %d (Press enter to terminate): " % count) # Inputting obsID's and counting up
		obsID_li.append(obsID) # Adding obsID's to the list
		count += 1
		if obsID == '':  # Enter breaks the loop
			del obsID_li[-1] # Removing the last addition to the list which is the enter --> ['####','####','']
			print obsID_li
			#obsID_checker()
			return obsID_li
			break

#quoted_list = obsID_selection()
#ordered_list = sorted(quoted_list, key=lambda x: float(x))
#print ordered_list

##################
### New Things ###
##################

def espace_filt(data, ccd):

	os.system("punlearn dmcopy")
	nrgy1 = raw_input("Input minimum energy filter value: ") # Allows specification of energy filter
	nrgy2 = raw_input("Input maximum energy filter value: ")
	for obs in data:
		os.system("dmcopy '%s_reproj_evt.fits[energy=%s:%s, ccd_id=%s]' %s_efilter.fits opt=all clobber=yes" % (obs, nrgy1, nrgy2, ccd, obs))

#espace_filt(ordered_list, ccd_id)

def extract_flare(data, time):

	for obs in data: # Cycles through obsID's
		os.system("dmlist %s_background.fits'[GTI7]' data" % obs)
		time_bin = []
		start = time[0]
		stop = time[1]
		bin = raw_input("Input bin length (usually 200): ")
		os.system("punlearn dmextract")
		os.system("dmextract infile=%s_background.fits'[bin time=%s:%s:%s]' outfile=%s_background.lc opt=ltc1 clobber=yes" % (obs, start, stop, bin, obs))
		os.system("deflare %s_background.lc %s_bkg_deflare.gti method=clean plot=yes save=%s_plot" % (obs, obs, obs))
	flare_checker() # Passes to the checker

#extract_flare(ordered_list, times)


# Finding the ccd_id from file header #

ccd_id = []
times = []

#def ccd_id_finder():
#
#	os.chdir(ENTER PATH HERE)
#	ccdlist = fits.open(ENTER FILE HERE broad_flux.fits)
#	ccd = ccdlist[0].header # Primary header
#	ID = ccd["CCD_ID"] # Get value of keyword 'CCD_ID', the ccd_id time of the observation
#	ccd_id.append(ID)
#	return ccd_id

# Finding the start and stop times from file header #

def time_finder():

	os.chdir("junk/reprojected_data")
	timelist = fits.open("4945_background.fits")
	tim = timelist[0].header # Primary header
	starttime = tim["TSTART"] # Get value of keyword 'TSTART', the start time of the observation
	stoptime = tim["TSTOP"] # Get value of keyword 'TSTOP', the start time of the observation
	times.extend([starttime, stoptime])
	return times

#ccd_id_finder()
#print ccd_id
time_finder()
print times

# TSTART and TSTOP do not seem to match the stop / start times from dmlist #
# CIAO does not have astropy #
