#!/usr/bin/env python

#######################################################
### Master code for the Summer project work of 2015 ###
#######################################################

 #####									    #####
##### By Dominic Eggerman and Grant Trembley #####
 #####									    #####

##########################################################
### Remember to launch CIAO before initiating the code ###
##########################################################


# I am replacing most of the CIAO code with print statements for now, later find and replace print with os.system.

 ###	   ###
### Imports ###
 ###	   ###

import numpy
import os
import sys
import glob
import pyfits
import subprocess
import os.path

###		###		###		###		###

 ###					  				  ###
### Directory creation and obsID selection ###
 ###					 				  ###

# Master directory creation #

def mkdir():
	directory = raw_input("Enter directory for work (ommit the last backslash '/') (e.g. ~/User/newobject): ")
	os.system("mkdir "+directory)
	os.path.expanduser(directory) # Help here ?

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

# Check for user satisfaction #

def obj_checker():

	choice = raw_input("Have you found the obsID's you wanted? (yes/no): ")
	for decision in choice:
		if choice == "yes":
			break
		elif choice == "no":
			find_obj() # Goes back to find_obj() function
			break
		else:
			choice = raw_input("Enter yes or no: ")

# Selecting obsID's #

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
			obsID_checker()
			return obsID_li
			break

# Converting, sorting and stringing obsID's #

def obsID_list_converter(obsIDs):

	obsID_list = map(int, obsIDs) # Converts all of the strings in the list to integers
	obsID_list.sort() # Sorting the list
	my_obsIDs = ",".join(map(str, obsID_list)) # Creates the string of obsID's for the download step
	return my_obsIDs

# Function to ask the user if the obsid selections are correct #

def obsID_checker():

	choice = raw_input("Are these the obsID's you wanted? (yes/no): ")
	for decision in choice:
		if choice == "yes":
			break
		elif choice == "no":
			obsID_li[:] = []
			obsID_selection()
			break
		else:
			choice = raw_input("Enter yes or no: ")

# Function activators 1 #

mkdir()
find_obj()
quoted_list = obsID_selection()
preped_list = obsID_list_converter(quoted_list)

###		###		###		###		###

 ###				   ###
### Downloading obsID's ###
 ###				   ###

# Download osbID's #

def chandra_download(string_list_to_download):

	print "download_chandra_obsid %s" % string_list_to_download

# Reduce data to lvl 2 #

def lvl2_data(lvl):

	for obs in lvl:
		print ("chandra_repro indir='{0}' outdir='{0}'/repro verbose=5 cleanup=no" .format(obs))

# Function activators 2 #

chandra_download(preped_list)
lvl2_data(quoted_list)

###		###		###		###		###

#########################################################################

# Aligning observations to the same WCS space and creating exposure corrected merged mosaic

def WCS():

	os.system("mkdir reprojected_data")
	print "reproject_obs '*/repro/*evt2.fits' reprojected_data/ verbose=5"
	os.system("mkdir exposure_corrected_mosaic")
	print "flux_obs 'reprojected_data/*reproj_evt.fits' exposure_corrected_mosaic/ bands=broad,csc bin=1 verbose=5"

# Opening broad_flux.img and waiting for the .reg files to be saved.

def deep_space_9():

	os.chdir("exposure_corrected_mosaic")
	print "ds9 broad_flux.img"  # Opens broad_flux.img in DS9 for making regions
	print "Create circular regions around the target and around a background portion of the image in DS9.  Save the target region as cluster.reg and the background region as background.reg.  Save both in /reprojected_data.  Make sure to save the regions in ciao format."
	raw_input("Press Enter to continue when the files have been saved...") # Allows a pause while the user creates the files
	os.chdir("../reprojected_data")

def file_checker(): # Checks for the .reg files, only runs twice as of now

	for fil in glob.glob("*.reg"):
		print "The files have been found"
		break
	else:
		print "Please create cluster.reg and background.reg in DS9"
		raw_input("Press Enter to continue when the files have been saved...")
		reg_files

# Function activators 2 #

file_checker()

###		###		###		###		###

###############################################################

# Sort quoted list, ordered_list is the new master list #

ordered_list = sorted(quoted_list, key=lambda x: float(x))

 ###																														  ###
### Filtering the reprojected files in energy space and extracting background only lightcurve by excluding cluster and sources ###
 ###																														  ###

# Filtering the files in energy space #

def espace_filt(data):

	print "punlearn dmcopy"
	nrgy = raw_input("Input energy filter values like so --> min:max : ") # Allows specification of energy filter
	ccd = raw_input("Input ccd_id: ") # Allows input of ccd_id
	for obs in data:
		print "dmcopy '%s_reproj_evt.fits[energy=%s, ccd_id=%s]' %s_efilter.fits opt=all clobber=yes" % (obs, nrgy, ccd, obs)

# Excluding cluster from the background #

def bkg_lightcurve(data):

	for obs in data:
		print "dmcopy '%s_efilter.fits[exclude sky=region(cluster.reg)]' %s_background.fits opt=all clobber=yes" % (obs, obs)

# Listing data to find the start and stop times #
# Needs to be automated

def extract_flare(data):

	for obs in data: # Cycles through obsID's
		print ("dmlist %s_background.fits'[GTI7]' data") % obs
		start = raw_input("Input desired START time from the table: ")
		stop = raw_input("Input desired STOP time from the table: ")
		bin = raw_input("Input bin length (usually 200): ")
		print ("punlearn dmextract")
		print ("dmextract infile=%s_background.fits'[bin time="+start+":"+stop+":"+bin+"]'"+" outfile=%s_background.lc opt=ltc1 clobber=yes") % (obs, obs)
		print ("deflare %s_background.lc %s_bkg_deflare.gti method=clean plot=yes save=%s_plot") % (obs, obs, obs)
	flare_checker() # Passes to the checker

# Checks with user to confirm the binning amount #

def flare_checker():
	choice = raw_input("Is the binning what you wanted? (yes/no): ")
	for decision in choice:
		if choice == "yes":
			break
		elif choice == "no": # Loops infinitly if you keep hitting "no" until a desired bin is found
			extract_flare(ordered_list)
			break
		else:
			choice = raw_input("Enter yes or no: ")

# Filter the event lists for all of the obsID's #

def evt_list_filt(data):

	for obs in data:
		print "dmcopy '%s_efilter.fits[@%s_bkg_deflare.gti]' %s_reproj_clean.fits clobber=yes" % (obs, obs, obs)
		# Left out the dmkeypar check exposure length step which was manual during unix work

# Function activators 3 #

espace_filt(ordered_list)
bkg_lightcurve(ordered_list)
extract_flare(ordered_list)
evt_list_filt(ordered_list)

###		###		###		###		###

 ###					 ###
### Blank Sky Backgrounds ###
 ###					 ###

# Create directories for the blank sky work #

def bsky_organiser():

	print "os.chdir('../'')"
	print "os.system('mkdir blank_sky')"
	print "cp *reproj_clean.fits "+directory+"/blank_sky/" # This will copy all of the ####_reproj_clean.fits files to the blank sky folder
	print "os.chdir('blank_sky/')"

# Get the background file #

def get_evt2():

	ID = ordered_list[0]
	print "loca = subprocess.check_output(['acis_bkgrnd_lookup', '%s_reproj_clean.fits'])" % ID
	# Above is the subprocess command, which captures the terminal output from the ciao background lookup and sets it = loca
	print "cp %s bkgevt2.fits" % ID # change ID for loca when putting in os commands ?
	print "dmcopy 'bkgevt2.fits[status=0]' bkgevt2_clean.fits" # Is status always = 0 ??
	# Not checking if GAINFILE is the same for background file and science files, is there a GAINFILE setter?

# Add pointing header keywords to the background file #

def evt2_pointer(data):

	for obs in data:
		print "dmmakepar %s_reproj_clean.fits %s_event_header.par" % (obs, obs)
		print "grep _pnt %s_event_header.par > %s_event_pnt.par" % (obs, obs)
		print "cp bkgevt2_clean.fits %s_bkgevt2_notproj.fits" % obs # Clone the clean background file into separate versions, one for each ObsID
	print "chmod +w *_bkgevt2_notproj.fits" # Make the clones writable
	for obs in data:
		print "dmreadpar %s_event_pnt.par '%s_bkgevt2_notproj.fits[events]' clobber+" % (obs, obs) # Migrate the pointing header keywords to the new clones

# While still in blank_sky, finding and copying aspect solution files over for reproject_events #

def aspect_sol(ID):

	for obs in ID:
		print "punlearn reproject_events"
		asp_file = [os.path.basename(x) for x in glob.glob('../%s/repro/*pcad*asol*')] % obs
		aspect = asp_file[0]
		print "cp ../%s/repro/*asol*.fits ." % obs
		print "reproject_events infile=%s_bkgevt2_notproj.fits outfile=%s_bkg_reproj_clean.fits aspect=%s match=%s_reproj_clean.fits random=0 verbose=5 clobber=yes" % (obs, obs, aspect, obs)
		del asp_file
		del aspect

# Function activators 4 #

bsky_organiser()
get_evt2()
evt2_pointer(ordered_list)
aspect_sol(ordered_list)

###		###		###		###		###

###############################################################################################

# Needs edits (mainly checks)

 ###			   ###
### Contour Binning ###
 ###			   ###

# Creating directory for work #

def contbin_dir():

	print "os.chdir('../')"
	print "mkdir contbin"
	print "os.chdir('contbin')"
	print "cp ../reprojected_data/merged_evt.fits ." # Copies merged-evt.fits over from reprojected_data
	print "ds9 merged_evt.fits" # Opens mergeed_evt.fits in ds9
	print "Make a box region around the target to obtain the x and y coordinates to input..."

min_vals = []

# Creating the region for contbinning work #

def reg_creation():

	xmin = raw_input("Input the minimum x coordinate (physical): ")
	xmax = raw_input("Input the maximum x coordinate (physical): ")
	ymin = raw_input("Input the minimum y coordinate (physical): ")
	ymax = raw_input("Input the maximum y coordinate (physical): ")
	min_vals.append(xmin)
	min_vals.append(ymin)
	return min_vals
	# Put something here to prevent bad inputs
	nrgy = raw_input("Input energy filter values like so --> min:max : ")
	# Maybe add raw_input checker
	print "dmcopy 'merged_evt.fits[energy=%s][bin x=%s:1, y=%s:1]' contbin_input.fits" % (nrgy, xminmax, yminmax)
	raw_input("Now save the region you created as contbin_mask.reg (ciao format).  Make sure to save it in the contbin folder in this projects parent directory (Press enter to continue)")
	# Need file checker before proceeding

# Farith step and preparation to make region files #

def farith():

	###### Need heasoft - maybe something to check that it is installed and need a way to call the program, also need it for contbin ######
	print "heasoft"
	print "farith contbin_input.fits 0 temp.fits MUL"
	print "farith temp.fits 1 allones.fits ADD"
	print "rm temp.fits"
	print "dmcopy 'allones.fits[sky=region(contbin_mask.reg)][opt full]' mask.fits"
	print "rm allones.fits"
	print "contbin --mask=mask.fits --sn=30 --smoothsn=3 --constrainfill --constrainval=2. contbin_input.fits"
	# May need user interaction here and a possible way to check for good data

# Creating .reg files #

def regions_sn30(values):

	print "mkdir regions_sn30"
	x = values[0] # Takes the elements from the min_val list to put in the final step
	y = values[1]
	print "make_region_files --minx=%s --miny=%s --bin=1 --outdir=regions_sn30 contbin_binmap.fits" % (x, y)

# Function activators 5 #

contbin_dir()
reg_creation()
farith()
regions_sn30(min_vals)

###		###		###		###		###

################################################################################################

 ###			 ###
### Spectral Maps ###
 ###			 ###

# Making the region list #

def region_list():

	print "os.chdir('regions_sn30/')"
	print "os.system('find . -name '*reg' > base.txt')"
	print "os.system('sed -i '' 's/\.reg//g' ~/SummerWork2015/Supertest/base.txt')" # This syntax for sed works with basic MAC sed (BSD)
	print "os.system('sed -i '' 's/\.\///g' ~/SummerWork2015/Supertest/base.txt')" # Removes the . / and .reg from the file names
	print "os.system('git download grants repository for combine.awk or paste it here')" # Downloads combine.awk ? This may need to go in BadDog repository
	print "os.system('awk -f combine.awk base.txt > region_list.txt')"

# Starting new window and calling heasoft and ciao #

def new_window(): # Need to start heasoft before ciao somehow in a new window or run it at the start?

	print "x"

# Create and prepare the spectral maps directory #

def spectral_create():

	print "os.chdir('../../')"
	print "os.system('mkdir spectral_maps')"
	print "os.chdir('spectral_maps')"
	print "os.system('cp ../contbin/regions_sn30/base.txt .')"

# Create the subfolders for spectral maps #

def data_create(ID):

	print "os.system('mkdir data')"
	print "os.system('mkdir data/prepared_spectral_data')"
	for obs in ID: # Copies the #### folders from the parent folder to data
		print "os.system('cp -r ../%s data/')" % obs
	print "os.system('cp ../blank_sky/*reproj* data/prepared_spectral_data')" # Copies the ####_reproj_clean.fits and ####_bkg_reproj_clean.fits

def regions_create(): # Prepares the regions folder for extract_contbin_spectra.py

	print "os.system('mkdir regions')"
	print "os.system('cp ../contbin/regions_sn30/combine.awk regions/')"
	print "os.system('cp ../contbin/contbin_mask.reg regions/')"
	print "os.system('cp ../contbin/regions_sn30/regions_list.txt regions/')"
	print "os.chdir('regions/')"
	print "os.system('find ../contbin/regions_sn30/ -name '*.reg' -exec cp {} . \;')" # Finds and copies all of the .reg files (uses find incase there are many)
	print "download Grant's code from github"

### Run extract_contbin_spectra.py

# Functions activators 6 #

region_list()
new_window()
spectral_create()
data_create(ordered_list)
regions_create()

###		###		###		###		###

###############
##### END #####
###############