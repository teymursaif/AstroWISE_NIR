
# VIRCAM photometry pipeline
# By Teymoor Saifollahi, Kapteyn Astronomical Institute, University of Groningen
# November 2017

# importing libraries

import sys, os
from astro.main.aweimports import *
from astro.config.startup import *
import pyfits
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

#module unload awe
#export PATH="/data/users/saifollahi/awe/Linux-redhat-7.3-x86_64/develop/astro/bin:$PATH"

def intro_msg() :
    print ('\n')
    print ('----------------------------------------------------------------------------------------------------')
    print ('VIRCAM Photometry Pipeline v1.0')
    print ('----------------------------------------------------------------------------------------------------')

def assign_input_params() :
    global instrument_name, filter_name, obs_start_date, obs_end_date
    instrument_name = 'VIRCAM'
    filter_name = 'ESO-Ks-0002'
    obs_start_date = [2012,1,1]
    obs_end_date = [2018,1,1]
    print ('- Input Parametes')
    print ('--- Instrument = ' + str(instrument_name))
    print ('--- Filter = ' + str(filter_name))
    print ('--- Starting Date of Observations = ' + str(obs_start_date))
    print ('--- Ending Date of Observations = ' + str(obs_end_date))

def retrieve_new_astrom(fits_file_name) :

    q = ((ReducedScienceFrame.instrument.name==instrument_name) & \
    (ReducedScienceFrame.filter.name==filter_name) & \
    (ReducedScienceFrame.creation_date>datetime.datetime(obs_start_date[0], \
    obs_start_date[1],obs_start_date[2])) & \
    (ReducedScienceFrame.creation_date<datetime.datetime(obs_end_date[0], \
    obs_end_date[1],obs_end_date[2]))).user_only()
    N = len(q)
    for i in range(0,N) :
        ap = (AstrometricParameters.reduced == q[i]).max('creation_date')
        if ap.RMS > 0.001:
            reduced = str(ap.reduced.filename)
            header = darma.header.header(ap.reduced.filename)
            header['CD1_1']=ap.CD1_1
            header['CD1_2']=ap.CD1_2
            header['CD2_1']=ap.CD2_1
            header['CD2_2']=ap.CD2_2
            header['CTYPE1']='RA---TAN'
            header['CTYPE2']='DEC--TAN'
            header['PV1_0']=ap.PV1_0
            header['PV1_1']=ap.PV1_1
            header['PV1_2']=ap.PV1_2
            header['PV1_3']=ap.PV1_3
            header['PV1_4']=ap.PV1_4
            header['PV1_5']=ap.PV1_5
            header['PV1_6']=ap.PV1_6
            header['PV2_0']=ap.PV2_0
            header['PV2_1']=ap.PV2_1
            header['PV2_2']=ap.PV2_2
            header['PV2_3']=ap.PV2_3
            header['PV2_4']=ap.PV2_4
            header['PV2_5']=ap.PV2_5
            header['PV2_6']=ap.PV2_6
            header.save(ap.reduced.filename+".head")
            print ('--- new astrometry has been saved on ' + ap.reduced.filename+".head")

def update_astrom_header_on_fits(fits_file_name) :
    head_file_name = fits_file_name + '.head'
    command = 'cp ' + str(fits_file_name) + ' image.fits'
    os.system(command)
    command = 'cp ' + str(head_file_name) + ' image.head'
    os.system(command)
    fits_file_name = 'image.fits'
    head_file_name = 'image.head'

    header = open(head_file_name)
    for lines in header:
        line = lines.split()
        i = 0
        for word in line :
            i = i+1
            if word=='CD1_1' :
                cd11 = line[i+1]
            if word=='CD1_2' :
                cd12 = line[i+1]
            if word=='CD2_1' :
                cd21 = line[i+1]
            if word=='CD2_2' :
                cd22 = line[i+1]
            if word=='PV1_0' :
                pv10 = line[i+1]
            if word=='PV1_1' :
                pv11 = line[i+1]
            if word=='PV1_2' :
                pv12 = line[i+1]
            if word=='PV1_3' :
                pv13 = line[i+1]
            if word=='PV1_4' :
                pv14 = line[i+1]
            if word=='PV1_5' :
                pv15 = line[i+1]
            if word=='PV1_6' :
                pv16 = line[i+1]
            if word=='PV2_0' :
                pv20 = line[i+1]
            if word=='PV2_1' :
                pv21 = line[i+1]
            if word=='PV2_2' :
                pv22 = line[i+1]
            if word=='PV2_3' :
                pv23 = line[i+1]
            if word=='PV2_4' :
                pv24 = line[i+1]
            if word=='PV2_5' :
                pv25 = line[i+1]
            if word=='PV2_6' :
                pv26 = line[i+1]
    header.close()

    fits = pyfits.open(fits_file_name)
    header = fits[0].header
    data = fits[0].data
    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC--TAN'
    header['CD1_1']= cd11
    header['CD1_2']= cd12
    header['CD2_1']= cd21
    header['CD2_2']= cd22
    header['PV1_0']= pv10
    header['PV1_1']= pv11
    header['PV1_2']= pv12
    header['PV1_3']= pv13
    header['PV1_4']= pv14
    header['PV1_5']= pv15
    header['PV1_6']= pv16
    header['PV2_0']= pv20
    header['PV2_1']= pv21
    header['PV2_2']= pv22
    header['PV2_3']= pv23
    header['PV2_4']= pv24
    header['PV2_5']= pv25
    header['PV2_6']= pv26
    fits.writeto('new_image.fits')
    command1 = 'rm ' + str(fits_name)
    command2 = 'mv new_image.fits ' + str(fits_name)
    os.system(command1)
    os.system(command2)

def make_list_of_reduced_science_frames() :
    return 'fits.list'

def photometry_main(list_of_reduced_science_frames) :
    fits_list = open(str(list_of_reduced_science_frames), 'r')
    for fits_file_name in fits_list :
        fits_file_name = fits_file_name.split()
        fits_file_name = fits_file_name[0]
        print ('\n- Photometry of frame ' + str(fits_file_name))
        print ('--- Retrieving reduced science frame')
        #retrieve_reduced_science_frames()
        print ('--- Retrieving weight frame')
        #retrieve_weight_frame()
        print ('--- Retrieving new astrometric parameters')
        retrieve_new_astrom(fits_file_name)
        print ('--- Updating astrometric parameters of fits file')
        #update_astrom_header_on_fits()
intro_msg()
assign_input_params()
fits_list = make_list_of_reduced_science_frames()
photometry_main(fits_list)
