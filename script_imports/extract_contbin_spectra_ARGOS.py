#!/usr/bin/env python

import numpy
import os
import sys
import glob
import pyfits

### Code to extract spectra from regions produced by contour binning ###

def make_spectrum(evt2name,region,specname):

    os.system('punlearn dmextract')
    os.system('dmextract %s"[sky=region(%s)][bin pi]" %s wmap="[energy=500:7000][bin det=8]" error=gaussian' % (evt2name,region,specname))

def make_background(bgevt2name,region,bgspecname):

    os.system('punlearn dmextract')
    os.system('dmextract %s"[sky=region(%s)][bin pi]" %s error=gaussian' % (bgevt2name, region, bgspecname))

def set_ardlib(bpixname):

    os.system('punlearn ardlib')
    os.system('acis_set_ardlib %s absolutepath=yes' % bpixname)

def getregioncounts(evt2name,regionfile):

    # Create image
    os.system('punlearn dmcopy')
    os.system('dmcopy "%s[sky=region(%s)][energy=500:7000][bin sky=1]" img.fits'% (evt2name,regionfile))

    os.system('punlearn dmstat')
    os.system('dmstat "img.fits[sky=region(%s)]" centroid=no' % regionfile)
    os.system('pget dmstat out_sum > temp.dat')
    ctsfile = open('temp.dat','r')
    counts = float(ctsfile.read())
    ctsfile.close()
    os.remove('temp.dat')
    os.remove('img.fits')
    
    return counts

def make_arf(specname,arfname,weightname,pbkfile,mskfile):

    os.system('punlearn mkwarf')
    os.system('mkwarf infile="%s[WMAP]" outfile=%s egridspec=0.3:9.5:0.005 threshold=0 weightfile=%s spectrumfile="" clobber=yes pbkfile="%s" dafile=CALDB mskfile="%s"' % (specname,arfname,weightname,pbkfile,mskfile))

def make_rmf(rmfout,spec,asol):

    os.system('punlearn mkacisrmf')
    os.system('mkacisrmf infile=CALDB outfile=%s energy=0.3:9.5:0.005 channel=1:1024:1 chantype=PI wmap="%s[WMAP]" gain=CALDB asolfile=%s' % (rmfout,spec,asol))

#def make_weightmap(evt2name,regionname,weightname):
    
#    os.system('punlearn dmcopy')
#    os.system('dmcopy "%s[sky=region(%s)][energy=500:7000][bin det=8]" %s' % (evt2name,regionname,weightname))

def add_spectra(spectra,specout):

    # Generate expr string
    exprstr = spectra[0]
    for spec in spectra[1:]:
        exprstr += '+%s' % spec

    # BACKSCAL keyword is taken from first spectrum
    hdulist = pyfits.open(spectra[0])
    backscalval = hdulist[1].header['BACKSCAL']
    hdulist.close()

    # Sum spectra
    os.system('punlearn mathpha')
    os.system('mathpha expr="%s" units="C" outfil="%s" exposure="CALC" areascal="NULL" backscal=%e ncomments=1 ERRMETH="POISS-0" comment1="Added files %s with MATHPHA"' % (exprstr,specout,backscalval,exprstr))

def add_arfs(arffiles,arfout,weights):

    # Generate expr string
    exprstr = arffiles[0]
    wgtstr = str(weights[0])
    for arf,wgt in zip(arffiles[1:],weights[1:]):
        exprstr += ',%s' % arf
        wgtstr += ',%f' % wgt
        
    os.system('punlearn addarf')
    os.system('addarf %s %s %s' % (exprstr,wgtstr,arfout))

def add_rmfs(rmffiles,rmfout,weights):

    # Generate expr string
    exprstr = rmffiles[0]
    wgtstr = str(weights[0])
    for rmf,wgt in zip(rmffiles[1:],weights[1:]):
        exprstr += ',%s' % rmf
        wgtstr += ',%f' % wgt
        
    os.system('punlearn addrmf')
    os.system('addrmf %s %s %s' % (exprstr,wgtstr,rmfout))

def generate_spectra(path,obsids,onchip,regionfile):

    weights = []
    weightmaps = []
    spectra = []
    bgspectra = []
    arffiles = []
    rmffiles = []
    sumweightmap = regionfile[:-4] + '_sum_detmap.fits'

    for obs in obsids:

        print "Processing obs. id %s ..." % obs
        obsidpath = path + str(obs)

        # Construct parameter list
        evtpath = path + 'prepared_spectral_data/'
        name = '%s_%s' % (regionfile[:-4],str(obs))

        params = {'evt2': evtpath + str(obs) +  '_reproj_clean.fits',
                  'bgevt2': evtpath + str(obs) + '_bkg_reproj_clean.fits',
                  'bpix': glob.glob(obsidpath + '/repro/acisf*_bpix1.fits')[0],
		  'reg': regionfile,
                  'spec': '%s_spec.pi' % name,
                  'bgspec': '%s_bgspec.pi' % name,
                  'rmf': '%s.rmf' % name,
                  'arf': '%s.arf' % name,
                  'weight': '%s.weight' % name,
                  'detmap': '%s_detmap.fits' % name,
                  'asol': glob.glob(obsidpath + '/repro/pcadf*_asol1.fits')[0],
                  'pbk': glob.glob(obsidpath + '/repro/acisf*_pbk0.fits')[0],
                  'msk': glob.glob(obsidpath + '/repro/acisf*_msk1.fits')[0],
                  }

        # Determine whether region is on chip, get region counts 0.5-7keV
        ctsweight = getregioncounts(params['evt2'],params['reg'])

        if ctsweight < 1.0:
            print "Region off chip, excluding this obs. id."
            continue
        
        # Make spectra!
        set_ardlib(params['bpix'])
        make_spectrum(params['evt2'],params['reg'],params['spec'])
        make_background(params['bgevt2'],params['reg'],params['bgspec'])
        make_arf(params['spec'],params['arf'],params['weight'],params['pbk'],params['msk'])
        make_rmf(params['rmf'],params['spec'],params['asol'])

        # Get weight for total arf
        weights.append(ctsweight)
        spectra.append(params['spec'])
        bgspectra.append(params['bgspec'])
        arffiles.append(params['arf'])
        rmffiles.append(params['rmf'])

    # Turn into arrays
    weights = numpy.array(weights)
    spectra = numpy.array(spectra)
    bgspectra = numpy.array(bgspectra)
    arffiles = numpy.array(arffiles)
    rmffiles = numpy.array(rmffiles)

    # Sum together those on same chip!
    #for i in numpy.array([0,1,2,3,7]):
    for i in numpy.array([7]):
   
        sclweights = weights[onchip == i]/sum(weights[onchip == i])

        print "Summing spectra ..."
        specout = regionfile[:-4] + '_sumc%d_spec.pi' % i
        add_spectra(spectra[onchip == i],specout)

        print "Summing background spectra ..."
        bgspecout = regionfile[:-4] + '_sumc%d_bgspec.pi' % i
        add_spectra(bgspectra[onchip == i],bgspecout)
    
        # Produce weighted arf file
        print "Weighting total arf by counts ..."
        arfout = regionfile[:-4] + '_sumc%d.arf' % i
        add_arfs(arffiles[onchip == i],arfout,sclweights)

        # Run mkacisrmf on total weight map
        print "Generating rmfs ..."
        rmfout = regionfile[:-4] + '_sumc%d.rmf' % i
        add_rmfs(rmffiles[onchip == i],rmfout,sclweights)

        # Group spectra
        print "Grouping spectrum ..."
        grpout = regionfile[:-4] + '_sumc%d_grp20.pi' % i
        os.system('punlearn grppha')
        os.system('grppha infile="%s" outfile="%s" chatter=0 comm="group min 20 & chkey BACKFILE %s & chkey RESPFILE %s & chkey ANCRFILE %s & exit"' % (specout,grpout,bgspecout,rmfout,arfout))

        # Update spec header
        os.system('punlearn dmhedit')
        os.system("""dmhedit %s filelist=none operation=add key=ANCRFILE value="'%s'" datatype=string""" % (specout,arfout))
        os.system('punlearn dmhedit')
        os.system("""dmhedit %s filelist=none operation=add key=RESPFILE value="'%s'" datatype=string""" % (specout,rmfout))
        os.system('punlearn dmhedit')
        os.system("""dmhedit %s filelist=none operation=add key=BACKFILE value="'%s'" datatype=string""" % (specout,bgspecout))

    # Remove files?
    for i in xrange(len(spectra)):
        os.remove(spectra[i])
        os.remove(bgspectra[i])
        os.remove(arffiles[i])
        os.remove(rmffiles[i])
        
    os.system('rm *.weight')

if __name__ == "__main__":

    for reg in regions:
        generate_spectra(path,obsids,onchip,reg)
        