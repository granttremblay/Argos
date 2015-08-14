#!/usr/bin/env python

#######################################################
### Master code for the Summer project work of 2015 ###
#######################################################

 #####									    #####
##### By Dominic Eggerman and Grant Trembley #####
 #####									    #####

######################################################################################
### Remember to launch Heasoft and CIAO (in that order) before initiating the code ###
######################################################################################

### Changed CIAO's python path to use Anaconda's python ###

 ###	   ###
### Imports ###
 ###	   ###

import numpy
import os
import sys
import glob
import subprocess
from astropy.io import fits
from script_imports import extract_contbin_spectra_ARGOS as extract
#import xspec

###		###		###		###		###

# Check if ciao, heasoft, and contbin are installed ??? #

 ###					  				  ###
### Directory creation and obsID selection ###
 ###					 				  ###

# Master directory creation #

directory = []

def mkdir():

	direct = raw_input("Enter directory for work (ommit the last backslash '/' and do not use ~, use /Users/usr/etc.) (e.g. /Users/usr/Project/target): ")
	make = os.system("mkdir "+direct)
	if os.path.exists(direct):
		os.chdir(direct)
		directory.append(direct)
	else:
		del direct
		del make
		mkdir()

# Loop for finding obsID's #

object_name = []

def find_obj():

	obj = raw_input("Enter object name (in quotes): ")
	os.system('find_chandra_obsID '+obj)
	new_obj = obj.replace("'", "").replace("\"","") # Replaces the '' or "" in the user input to check to see if ciao ran an error
	object_name.append(new_obj)
	ject = subprocess.Popen(['find_chandra_obsID', new_obj], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = ject.communicate()
	if ject.returncode != 0: # If an error, subprocess return code is not 0, this will pass again
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

# Sort quoted list, ordered_list is the new master list #

ordered_list = sorted(quoted_list, key=lambda x: float(x))

###		###		###		###		###

 ###				   ###
### Downloading obsID's ###
 ###				   ###

# Download osbID's #

def chandra_download(string_list_to_download):

	os.system('download_chandra_obsid '+string_list_to_download)

# Reduce data to lvl 2 #

def lvl2_data(lvl):

	for obs in lvl:
		os.system("chandra_repro indir=%s outdir=%s/repro verbose=5 cleanup=no" % (obs, obs))

# Function activators 2 #

chandra_download(preped_list)
lvl2_data(ordered_list)

###		###		###		###		###

 ###																						 ###
### Aligning observations to the same WCS space and creating exposure corrected merged mosaic ###
 ###																						 ###

# Align observations to the same WCS space.  Create directories and open broad_flux.fits #

def WCS():

	os.system('mkdir reprojected_data')
	os.system("reproject_obs '*/repro/*evt2.fits' reprojected_data/ verbose=5")
	os.system('mkdir exposure_corrected_mosaic')
	os.system("flux_obs 'reprojected_data/*reproj_evt.fits' exposure_corrected_mosaic/ bands=broad,csc bin=1 verbose=5")
	os.chdir('exposure_corrected_mosaic')
	os.system('cp broad_flux.img broad_flux.fits') # Copies broad_flux.img to broad_flux.fits so I can view it in my DS9
	print "Create circular regions around the target and around a background portion of the image in DS9.  Save the target region as cluster.reg and the background region as background.reg.  Save both in /reprojected_data/.  Make sure to save the regions in ciao format."
	raw_input("Press enter to open DS9, then close DS9 once you have fininshed making the files...")
	os.system('ds9 broad_flux.fits')  # Opens broad_flux.fits in DS9 for making regions
	os.chdir('../reprojected_data')

# Waiting for the .region files to be saved #

def deep_space_9():
	
	raw_input("Press Enter to continue when the files have been saved...") # Allows a pause while the user creates the files
	reg_file_check()
	
# Checks reprojected_data/ for the .region files #

def reg_file_check():

	if os.path.exists("./cluster.reg") == True:
		print "The file cluster.reg has been found"
	else:
		print "Please create cluster.reg in DS9 (save as ciao file) and save it to reprojected_data/"
		deep_space_9()
	if os.path.exists("./background.reg") == True:
		print "The file background.reg has been found"
	else:
		print "Please create background.reg in DS9 (save as ciao file) and save it to reprojected_data/"
		deep_space_9()

# Function activators 2 #

WCS()
deep_space_9()

###		###		###		###		###

 ###																														  ###
### Filtering the reprojected files in energy space and extracting background only lightcurve by excluding cluster and sources ###
 ###																														  ###

# Finding the ccd_id from file header #

ccd_id = []

def ccd_id_finder():

	os.chdir("../exposure_corrected_mosaic")
	ccdlist = fits.open("broad_flux.fits")
	ccd = ccdlist[0].header # Primary header
	ID = ccd["CCD_ID"] # Get value of keyword 'CCD_ID', the ccd_id time of the observation
	ccd_id.append(ID)
	os.chdir("../reprojected_data")
	return ccd_id

# Filtering the files in energy space #

def espace_filt(data, ccd):

	os.system("punlearn dmcopy")
	nrgy1 = raw_input("Input minimum energy filter value: ") # Allows specification of energy filter
	nrgy2 = raw_input("Input maximum energy filter value: ")
	choice = raw_input("Are these the values you wanted? (yes/no): ")
	for decision in choice:
		if choice == "yes":
			break
		elif choice == "no":
			del nrgy1
			del nrgy2
			espace_filt(ordered_list)
			break
		else:
			choice = raw_input("Enter yes or no: ")
	for obs in data:
		ccd_id = ccd[0]
		os.system("dmcopy '%s_reproj_evt.fits[energy=%s:%s, ccd_id=%s]' %s_efilter.fits opt=all clobber=yes" % (obs, nrgy1, nrgy2, ccd_id, obs))

# Excluding cluster from the background #

def bkg_lightcurve(data):

	for obs in data:
		os.system("dmcopy '%s_efilter.fits[exclude sky=region(cluster.reg)]' %s_background.fits opt=all clobber=yes" % (obs, obs))

# Finding the start and stop times from file header #

times = []

def time_finder(ccd):

	timelist = fits.open("4945_background.fits")
	ID = ccd[0]
	for num in range(2, 10):
		tim = timelist[num].header
		if tim["CCD_ID"] == ID:
			starttime = tim["TSTART"] # Get value of keyword 'TSTART', the start time of the observation ??
			stoptime = tim["TSTOP"] # Get value of keyword 'TSTOP', the start time of the observation
			times.extend([starttime, stoptime])
			return times
			break
		else:
			continue

# Deflaring #

def extract_flare(data, time):

	for obs in data: # Cycles through obsID's
		time_bin = []
		start = time[0]
		stop = time[1]
		bin = raw_input("Input bin length (usually 200): ")
		os.system("punlearn dmextract")
		os.system("dmextract infile=%s_background.fits'[bin time=%s:%s:%s]' outfile=%s_background.lc opt=ltc1 clobber=yes" % (obs, start, stop, bin, obs))
		os.system("deflare %s_background.lc %s_bkg_deflare.gti method=clean plot=yes save=%s_plot" % (obs, obs, obs))
	flare_checker() # Passes to the checker

# Checks with user to confirm the binning amount #

def flare_checker():
	choice = raw_input("Is the binning what you wanted? (yes/no): ")
	for decision in choice:
		if choice == "yes":
			break
		elif choice == "no": # Loops infinitly if you keep hitting "no" until a desired bin is found.
			os.system('rm *_bkg_deflare.gti')
			os.system('rm *_background.lc')
			extract_flare(ordered_list)
			break
		else:
			choice = raw_input("Enter yes or no: ")

# Filter the event lists for all of the obsID's #

def evt_list_filt(data):

	for obs in data:
		os.system("dmcopy '%s_efilter.fits[@%s_bkg_deflare.gti]' %s_reproj_clean.fits clobber=yes" % (obs, obs, obs))
		# Left out the dmkeypar check exposure length step which was manual during unix work

# Function activators 3 #

ccd_id_finder()
espace_filt(ordered_list, ccd_id)
bkg_lightcurve(ordered_list)
time_finder(ccd_id)
extract_flare(ordered_list, times)
evt_list_filt(ordered_list)

###		###		###		###		###

 ###					 ###
### Blank Sky Backgrounds ###
 ###					 ###

# Create directories for the blank sky work #

def bsky_organiser():

	os.chdir('../')
	os.system('mkdir blank_sky')
	os.system('cp reprojected_data/*reproj_clean.fits blank_sky') # This will copy all of the ####_reproj_clean.fits files to the blank sky folder
	os.chdir('blank_sky/')

# Get the background file #

def get_evt2():

	ID = ordered_list[0]
	loca = subprocess.check_output(['acis_bkgrnd_lookup', '%s_reproj_clean.fits' % ID]) # This is a subprocess command, which captures the terminal output from the ciao background lookup and sets it = loca
	location = loca[:-1] # There is a newline at the end of loca variable, have to remove it before passing to os.system
	os.system("cp %s bkgevt2.fits" % location)
	os.system("dmcopy 'bkgevt2.fits[status=0]' bkgevt2_clean.fits") # Is status always = 0 ??
	# Not checking if GAINFILE is the same for background file and science files, is there a GAINFILE setter?

# Add pointing header keywords to the background file #

def evt2_pointer(data):

	for obs in data:
		os.system("dmmakepar %s_reproj_clean.fits %s_event_header.par" % (obs, obs))
		os.system("grep _pnt %s_event_header.par > %s_event_pnt.par" % (obs, obs))
		os.system("cp bkgevt2_clean.fits %s_bkgevt2_notproj.fits" % obs) # Clone the clean background file into separate versions, one for each ObsID
	os.system("chmod +w *_bkgevt2_notproj.fits") # Make the clones writable
	for obs in data:
		os.system("dmreadpar %s_event_pnt.par '%s_bkgevt2_notproj.fits[events]' clobber+" % (obs, obs)) # Migrate the pointing header keywords to the new clones

# While still in blank_sky, finding and copying aspect solution files over for reproject_events #

def aspect_sol(ID):

	for obs in ID:
		os.system("punlearn reproject_events")
		asp_file = [os.path.basename(x) for x in glob.glob('../%s/repro/*pcad*asol*' % obs)] # Captures only the basename of the file in the path and adds it in a list
		aspect = asp_file[0] # Assigns the file to a variable
		os.system("cp ../%s/repro/*asol*.fits ." % obs)
		os.system("reproject_events infile=%s_bkgevt2_notproj.fits outfile=%s_bkg_reproj_clean.fits aspect=%s match=%s_reproj_clean.fits random=0 verbose=5 clobber=yes" % (obs, obs, aspect, obs))
		del asp_file
		del aspect

# Function activators 4 #

bsky_organiser()
get_evt2()
evt2_pointer(ordered_list)
aspect_sol(ordered_list)

###		###		###		###		###

 ###			   ###
### Contour Binning ###
 ###			   ###

# Creating directory for work #

def contbin_dir():

	os.chdir('../')
	os.system("mkdir contbin")
	os.chdir('contbin')
	os.system("cp ../reprojected_data/merged_evt.fits .") # Copies merged-evt.fits over from reprojected_data to the current folder
	print "Make a box region around the target to obtain the x and y coordinates to input (physical coordinates in ds9)..."
	print "For now, write these coordinates down. Make sure to save the region as contbin_mask.reg (ciao format) in the contbin folder in this project's parent directory"
	os.system("ds9 merged_evt.fits") # Opens merged_evt.fits in ds9
	
# Define coordinate boundaries #

min_values = []
energy_li = []

def coordinate_inp():

	xmin = raw_input("Input the minimum x coordinate (physical): ")
	xmax = raw_input("Input the maximum x coordinate (physical): ")
	ymin = raw_input("Input the minimum y coordinate (physical): ")
	ymax = raw_input("Input the maximum y coordinate (physical): ")
	min_values.extend([xmin, xmax, ymin, ymax]) # Adds all of the entered values to the list
	minmax_check(min_values)
	return min_values

# Check the min / max values #

def minmax_check(value): # Ask if the values are good for the user??

	for val in value:
		if value[0] < value[1]: # Checks if the max / min values have been correctly entered (that the max is not smaller than the min)
			pass
		else:
			print "The min / max values were entered incorrectly, please enter them again"
			del min_values[:]
			coordinate_inp()
			break
		if value[2] < value[3]:
			pass
		else:
			print "The min / max physical coordinate values were entered incorrectly, please enter them again"
			del min_values[:]
			coordinate_inp()
			break

# Define energy boundaries #

def energy_inp():

	nrgy1 = raw_input("Input minimum energy filter value: ")
	nrgy2 = raw_input("Input maximum energy filter value: ")
	energy_li.extend([nrgy1, nrgy2])
	nrgy_check(energy_li)
	return energy_li

# Checks the energy input values #

def nrgy_check(value): # Ask if the values are good for the user ??

	for val in value:
		if value[0] < value[1]:
			pass
		else:
			print "The min / max energy filter values were entered incorrectly, please enter them again"
			del energy_li[:]
			energy_inp()
			break

# Creating the region for contbinning work #

def reg_creator(energy, values): # DS9 closes before this step ??

	raw = energy + values
	minmaxnrgy_check(energy_li, min_values)
	os.system("dmcopy 'merged_evt.fits[energy=%s:%s][bin x=%s:%s:1, y=%s:%s:1]' contbin_input.fits clobber=yes" % tuple(raw))

# Asks if inputs are good #

def minmaxnrgy_check(energy, values):

	raw = energy + values
	print "min_energy:%s max_energy:%s min_x:%s max_x:%s min_y:%s max_y:%s" % tuple(raw)
	choice = raw_input("Are these the values you wanted? (yes/no): ")
	for decision in choice:
		if choice == "yes":
			break
		elif choice == "no":
			del min_values[:]
			del energy_li[:]
			coordinate_inp() # Goes back to coordinate_inp() function
			energy_inp()
			reg_creator(energy_li, min_values)
			break
		else:
			choice = raw_input("Enter yes or no: ")

# Checks for contbin_mask.reg #

def contbinmask_file_check(): # Checks for the .reg files, only runs twice as of now

	if os.path.exists("contbin_mask.reg") == True:
		print "The file has been found"
	else:
		print "Please create contbin_mask.reg in DS9 (save as ciao file to the contbin folder in the project directory)"
		raw_input("Press enter once the file has been made...")
		contbinmask_file_check()

# Farith step and preparation to make region files #

def farith():

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

# Making the region lists in base.txt #

def region_list(direct):

	os.chdir('regions_sn30')
	os.system("find . -name '*reg' > base.txt")
	for path in direct:
		os.system("sed -i '' 's/\.reg//g' %s/contbin/regions_sn30/base.txt" % path) # This syntax for sed works with basic MAC sed (BSD)
		os.system("sed -i '' 's/\.\///g' %s/contbin/regions_sn30/base.txt" % path) # Removes the . / and .reg from the file names

# Function activators 5 #

contbin_dir()
coordinate_inp()
energy_inp()
reg_creator(energy_li, min_values)
contbinmask_file_check()
farith()
regions_sn30(min_values)
region_list(directory)

###		###		###		###		###

################################################################################################

 ###			 ###
### Spectral Maps ###
 ###			 ###

# Needs edits before running ?? #

# Create and prepare the spectral maps directory #

def spectral_create():

	os.chdir('../../')
	os.system('mkdir spectral_maps')
	os.chdir('spectral_maps')
	os.system('cp ../contbin/regions_sn30/base.txt .')

# Create the subfolders for spectral maps #

def data_create(ID):

	os.system('mkdir data')
	os.system('mkdir data/prepared_spectral_data')
	for obs in ID: # Copies the obsID folders from the parent folder to data/
		os.system('cp -r ../%s ./data' % obs)
	os.system('cp ../blank_sky/*reproj* ./data/prepared_spectral_data') # Copies the ####_reproj_clean.fits and ####_bkg_reproj_clean.fits to the prepared_spectral_data folder.

# Prepares the regions folder for extract_contbin_spectra.py #

def regions_create(directory):

	os.system('mkdir regions')
	reg_list_write() # Passes to the file writer below
	os.system("cp ../../contbin/contbin_mask.reg .")
	for path in directory:
		os.system("find %s/contbin/regions_sn30 -name '*.reg' -exec cp {} . \;" % path) # Finds and copies all of the xaf_###.reg files to regions/ (uses find incase there are many)

# Writes the region_list.txt file #

def reg_list_write():

	with open ("base.txt", "r") as names:
		int_list = names.readlines() # Gets all of the xaf_###/n names from base.txt and assigns them to a list
	xaf = []
	for item in int_list:
		xaf_item = item[:-1] # Gets rid of the /n from the xaf_### names
		xaf.append(xaf_item)
	os.chdir("regions")
	combine = open("region_list.txt", "a") # Makes the file region_list.txt
	for x in xaf:
		combine.write("%s %s_sumc7_spec.pi\n" % (x, x)) # Writes the file with the xaf_### names and the ending that was previously done with combine.awk
	combine.close()

# Run extract_contbin_spectra.py #

def extract_spectra(direct, obs, ccd):

	# Heasoft should be launched before ciao
	# Set "path" equal to the user specified directory for work plus the /spectral_maps/data/ folder
	path = direct[0] + '/spectral_maps/data/'
	# ObsID's get passed as a list
	# Convert the ccd_id value to an intger and make a list as long as the obsID list ?? There may need to be different ccd_id's
	onchip = []
	chip_single = map(int, ccd)
	single = chip_single[0]
	for ID in obs:
		onchip.append(single)
	# Get the list of regions
	regions = glob.glob('xaf_*.reg') # Get list of regions
	# Call extract_contbin_spectra.py for each region. extract_contbin_spectra is set as alias extract
	for reg in regions:
		extract.generate_spectra(path,obs,onchip,reg)

# Functions activators 6 #

spectral_create()
data_create(ordered_list)
regions_create(directory)
extract_spectra(directory, ordered_list, ccd_id)

###		###		###		###		###

###############################################################################################

 ###	 ###
### Xspec ###
 ###	 ###

### Note: PyXspec works via "import xspec" on Linux but not on Mac for now ###

 # Getting data from NASA Extragalactic Database #

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

 # Pass the NED data to Xspec #

 # Make .tcl script in python

redshift_finder(object_name)

###############
##### END #####
###############