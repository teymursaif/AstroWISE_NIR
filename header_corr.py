
import sys, os
import numpy as np
import pyfits

# this code use a dark frame to produce a fake bias frame.
fits_name = sys.argv[1]
fits = pyfits.open(fits_name)
header = fits[0].header
data = fits[0].data
#keyword['PV2_1'] = 1.0
#keyword['PV2_2'] = 0.0
#keyword['PV2_3'] = 44.0
#keyword['PV2_4'] = 0.0
#keyword['PV2_5'] = -10300.0
header['CTYPE1'] = 'RA---TAN'
header['CTYPE2'] = 'DEC--TAN'
#header['CDELT1'] = 0.00009444
#header['CDELT2'] = 0.00009444
fits.writeto('newimage.fits')
command1 = 'rm ' + str(fits_name)
command2 = 'mv newimage.fits ' + str(fits_name)
os.system(command1)
os.system(command2)
fits.close()
