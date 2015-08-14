#!/usr/bin/env python
# File for testing code segments

import numpy
import os
import sys
import glob
import subprocess
from astropy.io import fits
from script_imports import extract_contbin_spectra_ARGOS as extract
from astroquery.ned import Ned

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

object_name = []

def find_obj():

	obj = raw_input("Enter object name (in quotes): ")
	new_obj = obj.replace("'", "").replace("\"","") # Replaces the '' or "" in the user input to check to see if ciao ran an error
	object_name.append(new_obj)

find_obj()

##################
### New Things ###
##################

NED_data = []

def redshift_finder(objname):

	obj = objname[0]
	main_table = Ned.query_object(obj)
	redshift = main_table['Redshift'][0]
	NED_data.append(redshift)
	RA = main_table['RA(deg)'][0]
	NED_data.append(RA)
	DEC = main_table['DEC(deg)'][0]
	NED_data.append(DEC)

redshift_finder(object_name)
