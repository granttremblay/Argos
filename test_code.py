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

quoted_list = obsID_selection()
ordered_list = sorted(quoted_list, key=lambda x: float(x))
print ordered_list

##################
### New Things ###
##################

min_values = []

def coordinate_inp():

	xmin = raw_input("Input the minimum x coordinate (physical): ")
	xmax = raw_input("Input the maximum x coordinate (physical): ")
	ymin = raw_input("Input the minimum y coordinate (physical): ")
	ymax = raw_input("Input the maximum y coordinate (physical): ")
	min_values.extend([xmin, xmax, ymin, ymax]) # Adds all of the entered values to the list
	return min_values

def farith():

	os.chdir("junkdir/contbin")
	###### Need heasoft - maybe something to check that it is installed and need a way to call the program, also need it for contbin ######
	os.system("farith contbin_input.fits 0 temp.fits MUL")
	os.system("farith temp.fits 1 allones.fits ADD")
	os.system("rm temp.fits")
	os.system("dmcopy 'allones.fits[sky=region(contbin_mask.reg)][opt full]' mask.fits")
	os.system("rm allones.fits")
	os.system("contbin --mask=mask.fits --sn=30 --smoothsn=3 --constrainfill --constrainval=2. contbin_input.fits")
	# May need user interaction here and a possible way to check for good data

# Creating .reg files #

def regions_sn30(values):

	os.system("mkdir regions_sn30")
	x = min_values[0] # Takes the elements from the min_val list to put in the final step
	y = min_values[2]
	os.system("make_region_files --minx=%s --miny=%s --bin=1 --outdir=regions_sn30 contbin_binmap.fits" % (x, y))

# Function activators 5 #

coordinate_inp()
farith()
regions_sn30(min_values)
