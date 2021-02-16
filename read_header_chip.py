# march 27th, 2017

import sys, os
import numpy as np
import pyfits 
import shutil

# getting the list of fits files in relevant directory.
fits_list = open('fits_list+.cat', 'r')

# defining desired keywords to read from header
keywords = ['DATE', 'EXPTIME', 'HIERARCH ESO DET CHIP ID', \
           'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CRPIX1', 'CRPIX2']

# defining more keywords

#
print "\n"
print keywords,
print "\n"

i = 0
# reading the list of fits files
for fits_filename in fits_list :

  fits_filename = fits_filename.split()
  fits_filename = fits_filename[0]
  fits_filename = "VIRCAM-DATA/" + fits_filename
  fits_file = pyfits.open(str(fits_filename)) 
  
  for i in range(1,17) :
  
    print "CHIP" + str(i) + "	",
  
    for keyword in keywords:

      data = fits_file[i].header[keyword]
      print fits_file[i].header[keyword], "	",

    print ''
  print ''
  
  fits_file.close()
  


