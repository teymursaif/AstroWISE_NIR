import sys, os
import numpy as np

fits_list = open('fits.list', 'r')
swarp_list = open('swarp.list', 'w')
for fits_file in fits_list :
    fits_id = fits_file[63:73]
    #print fits_id
    weight_list = open('weight.list', 'r')
    for weight_file in weight_list :
        weight_id = weight_file[63:73]
        if weight_id == fits_id :
            #print (fits_id, weight_id)
            weight_file = weight_file.split()
            fits_file = fits_file.split()
            command = 'cp ' + str(weight_file[0]) + '.fits '+ str(fits_file[0]) + '.weight.fits'
            print command
            #os.system(command)
            swarp_list.write(str(fits_file[0]) + '.fits' + '\n')
	    command2 = 'cp ' + str(fits_file[0]) + '.fits.head '+ str(fits_file[0]) + '.head'
            os.system(command2)
    weight_list.close()
fits_list.close()
swarp_list.close()
