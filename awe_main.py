import sys, os
from astro.main.aweimports import *
from astro.config.startup import *
import pyfits

#command = 'awe make_list.py'
#os.system(command)

#raw_list = open('raw_list.cat', 'r')
i = 0
raw_list = open(str(sys.argv[1]), 'r')
for raw in raw_list :
    i = i + 1
    if i> 826 and i < 1000 :
        raw=raw.split()
        raw=raw[0]
        sample=list()
        sample.append(raw)
        command = 'awe red.py ' + str(raw)
        os.system(command)
    """
    back = (BackgroundFrame.raw.filename==sample)
    if back :
        print ("BackgroundFrame does exist in DB.")
    else :
        print ("BackgroundFrame does not exist in DB. Creating BackgroundFrame ...")
        lpu.run('Background', i='VIRCAM', raw_filenames=sample, C=1)
    lpu.run('Reduce', i='VIRCAM', filter='ESO-Ks-0002', raw_filenames=sample, oc=0, C=1)
    os.system('rm tmp*')
    """
raw_list.close()

"""
cra,cdec = 54.5 , -35.5

raw = (str(sys.argv[1]))
sample = list()
sample.append(raw)
print (raw)
frames=(ReducedScienceFrame.filename==raw)

#back = (BackgroundFrame.raw.filename==raw)
#print (len(back))
#if back :
#    print ("BackgroundFrame does exist in DB.")
#else :
#    print ("BackgroundFrame does not exist in DB. Creating BackgroundFrame ...")
#    lpu.run('Background', i='VIRCAM', raw_filenames=sample, C=1)

#lpu.run('Reduce', i='VIRCAM', raw_filenames=sample, C=1)

for frame in frames :
    print (frame.filename)
    rmcommand = "rm " + frame.filename
    os.system(rmcommand)
    frame.retrieve()
    fitsname = "/data/users/saifollahi/awe/develop/"+str(frame.filename)
    framename = frame.filename

fits_file = pyfits.open(fitsname, mode='update')

keyword = fits_file[0].header

keyword['PV2_1'] = 1.0
keyword['PV2_2'] = 0.0
keyword['PV2_3'] = 44.0
keyword['PV2_4'] = 0.0
keyword['PV2_5'] = -10300.0
newimage = "C-"+framename
fits_file.writeto(newimage)
#fits_file.flush()
fits_file.close()
mvcommand = "mv " + newimage + " " + framename
os.system(rmcommand)
os.system(mvcommand)

lpu.run('Astrometry', i='VIRCAM', red_filenames=sample, C=0)

#ra,cdec = 54.5 , -35.5
#lpu.run('Regrid', i='VIRCAM', red_filenames=sample, gra=cra, gdec=cdec, C=0)

#q = ((RegriddedFrame.reduced.filename==raw)).max('creation_date')
#print (len(q))
#for i in q :
#    i.retrieve()
"""
