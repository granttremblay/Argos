#!/usr/bin/env python
# File for testing code segments

import numpy
import os
import sys
import glob
import pyfits
import subprocess
from astropy.io import fits
from script_imports import extract_contbin_spectra_ARGOS as extract
import xspec

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

##################
### New Things ###
##################

print "imported"