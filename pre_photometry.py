import sys, os
#from astro.main.aweimports import *
#from astro.config.startup import *
import pyfits
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.stats import norm
import scipy.interpolate

command = 'export PATH="/data/users/saifollahi/awe/Linux-redhat-7.3-x86_64/develop/astro/bin:$PATH"'
os.system(command)
list_name = 'fits.list'

chips=[22]
#chips = [35,36,38,39,41,42,43,44,45,46,47,22,23,25,30,33]
for chip in chips:
    chip='Virgo'+str(chip)
    chip_fits_file = open(str(chip)+'.list','w')
    directory = ''
    fits_list = open(list_name, 'r')
    for fits_file_name in fits_list :
        fits_file_name = fits_file_name.split()
        fits_file_name = fits_file_name[0]

        #
        head_file_name = fits_file_name + '.head'
        header = open(head_file_name)
        for lines in header:
            line = lines.split()
            i = 0
            for word in line :
                i = i+1
                if word=='DATE' :
                    date = line[i+1]
                    year = int(date[1:5])
                    month = int(date[6:8])
                    day = int(date[9:11])

                    #print year, month
        header.close()
        #
        if chip in fits_file_name :
        #and ( (year == 2013 and month == 1 and day == 30 ) ) :
            #print (year, month, day)
            #print fits_file_name
            #command1 = 'cp ../data+/'+str(fits_file_name) + ' .'
            #command2 = 'cp ../data+/'+str(fits_file_name)+'.weight.fits' + ' .'
            #command3 = 'cp ../data+/'+str(fits_file_name) +'.head' + ' .'
            #os.system(command1)
            #os.system(command2)
            #os.system(command3)
            chip_fits_file.write(fits_file_name+'\n')
    chip_fits_file.close()
    fits_list.close()

for chip in chips:
    chip='Virgo'+str(chip)
    #command1 = 'python photometry.py ' + str(chip)+ '.list' + ' ' + str(chip)
    command1 = 'python frame_photometry.py ' + str(chip)+ '.list' + ' ' + str(chip)
    #command2 = 'rm Sci*' + str(chip) + '*'
    print ('----------------' + str(chip) + '----------------\n')
    #print (chip)
    #print command

    os.system(command1)
    #os.system(command2)
