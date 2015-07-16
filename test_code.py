#!/usr/bin/env python
# File for testing code segments

import numpy
import os
import sys
import glob
import pyfits
import subprocess
import os.path

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
			# check = obsID_checker()
			return obsID_li
			break


#quoted_list = obsID_selection()
#ordered_list = sorted(quoted_list, key=lambda x: float(x))


### New things below here ###

# Loop for finding obsID's #

def find_obj():

	obj = raw_input("Enter object name (in quotes): ")
	os.system('find_chandra_obsID '+obj)
	new_obj = obj.replace("'", "").replace("\"","") # Replaces the '' or "" in the user input to check to see if ciao ran an error
	ject = subprocess.Popen(['find_chandra_obsID', new_obj], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = ject.communicate()
	if ject.returncode != 0: # If an error, subprocess return code is not 0, this will pass to redo()
		print "The object name you entered is not in the CIAO database. Enter the object ID again (in quotes)"
		find_obj()
	else:
		pass
	obj_checker()

def obj_checker():

	choice = raw_input("Have you found the obsID's you wanted? (yes/no): ")
	for decision in choice:
		if choice == "yes":
			break
		elif choice == "no":
			find_obj()
			break
		else:
			choice = raw_input("Enter yes or no: ")

find_obj()
