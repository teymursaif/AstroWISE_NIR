# march 27th, 2017

import sys, os
import numpy as np
import pyfits
import shutil

# getting the list of fits files in relevant directory.
#os.system("ls VIRCAM-DATA-2013/ > fits_list_2013.cat")
fits_list = open('mef.cat', 'r')

# defining desired keywords to read from header
keywords = ['ARCFILE', 'DATE', 'TELESCOP', 'INSTRUME', 'OBJECT', 'EXPTIME']

# defining more keywords

#keywords.append('HIERARCH ESO DPR TYPE')
#keywords.append('HIERARCH ESO INS FILT1 DATE')
#keywords.append('HIERARCH ESO INS FILT1 FOCUS')
#keywords.append('HIERARCH ESO INS FILT1 ID')
#keywords.append('HIERARCH ESO INS FILT1 NAME')
#keywords.append('HIERARCH ESO INS FILT1 NO')


# defining the output file (header_table)
header_table = open('header_table.cat','w')
os.system(" > header_table.cat")

header_table.write(str(keywords))
header_table.write("\n")
print "\n"
print keywords,
print "\n"

i = 0
# reading the list of fits files
for fits_filename in fits_list :
  i = i + 1
  print i,
  fits_filename = fits_filename.split()
  fits_filename = fits_filename[0]
  fits_filename = "VIRCAM-DATA-Ks/" + fits_filename
  fits_file = pyfits.open(str(fits_filename))

  for keyword in keywords:

    data = fits_file[0].header[keyword]
    header_table.write(str(data)+"	")
    print fits_file[0].header[keyword], "	",

  print ''

  #print str(fits_file[0].header.keys())
  #print "\n"

  header_table.write("\n")
  fits_file.close()

header_table.close()

# creating directories based on data_types
header_table = open('header_table.cat', 'r')
data_types = ['DARK', 'FLAT,TWILIGHT', 'STD,FLUX', 'Fornax', 'DARK,GAIN', 'FLAT,LAMP,GAIN', \
              'DARK,LINEARITY', 'DARK,CHECK', 'FLAT,LAMP,LINEARITY', 'FLAT,LAMP,CHECK' , 'NGC1399', 'NGC1399_K_t1_v01_nx']

for data_type in data_types :

  if not os.path.exists(data_type) :
    os.makedirs(data_type)

print "\nCopying files ...\n "
for lines in header_table :

  line = lines.split()

  if line[4] in data_types :

    data_type = str(line[4])
    data_arcfile = str(line[0])
    src = "VIRCAM-DATA-Ks"+data_arcfile
    des = str(data_type)+"/"
    #print src + " to " + des
    shutil.copy(src,des)
    print src

# creating PURPOSE directories + copying
if not os.path.exists("TWILIGHT") :
  os.makedirs("TWILIGHT")

if not os.path.exists("PHOTOM") :
  os.makedirs("PHOTOM")

os.system("cp FLAT,TWILIGHT/* TWILIGHT/")
os.system("cp NGC1399_K_t1_v01_nx/* PHOTOM/")
os.system("cp NGC1399/* PHOTOM/")
#os.system("cp STD,FLUX/* PHOTOM/")

#env project='VIRCAM@VISTA (4)' AWEVERSION='current' awe /Software/users/astro-wise/awe/current/astro/toolbox/ingest/Ingest.py -i *.fits -t photom -commit
#env project='VIRCAM@VISTA (4)' AWEVERSION='current' awe /Software/users/astro-wise/awe/current/astro/toolbox/ingest/Ingest.py -i *.fits -t twilight -commit
