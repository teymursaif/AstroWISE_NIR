#############################################################
######################## V I R C A M ########################
###################### P I P E L I N E ######################
#############################################################

First import :
import sys, os
from astro.main.aweimports import *
from astro.config.startup import *

Selecting the data
f = (RawTwilightFlatFrame.instrument.name == 'VIRCAM') & \
    (RawTwilightFlatFrame.DATE_OBS < datetime.datetime(2013,1,1)) & \
    (RawTwilightFlatFrame.DATE_OBS > datetime.datetime(2012,1,1)) & \

     do not forget about filter !

Inspecting the header and statistics of selected data :
for r in f: dir(r)
['DATE', 'DATE_OBS', 'EXPTIME', 'LST', 'MEAN_HIGH', 'MEAN_LOW', 'MJD_OBS', 'NAXIS1', 'NAXIS2', 'OBJECT', 'OBSERVER', 'OVSCX', 'OVSCXPRE', 'OVSCXPST', ...]

printing statistics of data
for r in f: print(r.DATE_OBS, r.filter.name, r.imstat.max, r.imstat.min, r.imstat.mean, r.imstat.stdev)

#############################################################
##################### R E D U C T I O N #####################
#############################################################

0. Cold pixel map :
lpu.run('ColdPixels', instrument='VIRCAM', filter='ESO-Ks-0002', date='2012-01-01', commit=False)

1. Background estimation :
lpu.run('Background', i='VIRCAM', raw_filenames=sample, C=0)

2. Creating twilight flat
lpu.run('TwilightFlat', instrument='VIRCAM', filter='ESO-Ks-0002', date='2010-01-01', commit=False)

3. Creating master flat from twilight flats
lpu.run('MasterFlat', instrument='VIRCAM', filter='ESO-Ks-0002', date='2012-01-01', combine=3, commit=False)

4. Reduction

f = (RawScienceFrame.instrument.name == 'VIRCAM') & \
    (RawScienceFrame.DATE_OBS < datetime.datetime(2018,11,1)) & \
    (RawScienceFrame.DATE_OBS > datetime.datetime(2012,1,1) &\
    (RawScienceFrame.filter.name == 'ESO-J-0002')).user_only()

import sys, os
for i in f :
    sample=list()
    sample.append(i.filename)
    back = (BackgroundFrame.raw.filename==i.filename)
    if back :
        print ("BackgroundFrame does exist in DB.")
    else :
        print ("BackgroundFrame does not exist in DB. Creating BackgroundFrame ...")
        lpu.run('Background', i='VIRCAM', raw_filenames=sample, C=1)
    lpu.run('Reduce', i='VIRCAM', filter='ESO-Ks-0002', raw_filenames=sample, oc=0, C=1)
    os.system('rm tmp*')

#############################################################
#################### A S T R O M E T R Y ####################
#############################################################
# first ZPN -> TAN in header of frames and then performing Astrometry
5. local Astrometry
First, optimize sextractor parameters ! -> config.py

module unload awe
export PATH="/data/users/saifollahi/awe/Linux-redhat-7.3-x86_64/develop/astro/bin:$PATH"
q = ((ReducedScienceFrame.instrument.name=='VIRCAM') & \
(ReducedScienceFrame.filter.name=='ESO-Ks-0002') &\
(ReducedScienceFrame.creation_date>datetime.datetime(2018,3,1)))
sample = list()
for i in q :
    sample = list()
    sample.append(i.filename)
    #i.retrieve()
    command='python header_corr.py '+ str(i.filename)
    os.system(command)
    lpu.run('Astrometry', i='VIRCAM', red_filenames=sample, C=1)
    command='rm tmp*'
    os.system(command)
    command='rm *Virgo*'
    os.system(command)

#Astrometry Test
import os, sys
fits_name = 'Sci-TSAIFOLLAHI-VIRCAM-------ESO-Ks-0002-ESO-Virgo33-Red---Sci-58031.6212486-a216f9561ccca7ce5ffe3cb9646fcb2393f4c467.fits'
command='python header_corr.py '+ str(fits_name)
os.system(command)
sample = ['Sci-TSAIFOLLAHI-VIRCAM-------ESO-Ks-0002-ESO-Virgo33-Red---Sci-58031.6212486-a216f9561ccca7ce5ffe3cb9646fcb2393f4c467.fits']
dpu.run('Astrometry', i='VIRCAM', red_filenames=sample, C=0)
#i.retrieve()

6. Creating global source list

q = ((AstrometricParameters.instrument.name=='VIRCAM') & \
(AstrometricParameters.quality_flags == 0) & \
(AstrometricParameters.is_valid >= 1) & \
(AstrometricParameters.filter.name=='ESO-Ks-0002') & \
(AstrometricParameters.creation_date > datetime.datetime(2017,5,14))).user_only()
sample=[]
for i in q :
    sample.append(i.reduced.filename)

task = GAstromSourceListTask(red_filenames=sample, commit=0)
task.execute()

7. global Astrometry (SKIPPING!)
task = GAstromTask(red_filenames=sample, C=0)
task.execute()

#############################################################
#################### P H O T O M E T R Y ####################
#############################################################
8. Photometry

q = ((AstrometricParameters.instrument.name=='VIRCAM') & \
(AstrometricParameters.filter.name=='ESO-Ks-0002') & \
(AstrometricParameters.quality_flags == 0) & \
(AstrometricParameters.is_valid >= 1) & \
(AstrometricParameters.creation_date > datetime.datetime(2017,5,14))).user_only()
sample=[]
for i in q :
    sample.append(i.reduced.filename)

# creating source catalogues
from astro.recipes.PhotCalExtractResulttable import PhotcatTask
task = PhotcatTask(instrument = 'VIRCAM', red_filenames = sample, C=0)
task.execute()

from astro.recipes.PhotCalExtractZeropoint import PhotomTask
task = PhotomTask(instrument = 'VIRCAM', red_filenames = sample, C=0)
task.execute()

#############################################################
#################### C O A D D I T I O N ####################
#############################################################
9. Regrid

cra,cdec = 54.5 , -35.5
q2 = ((AstrometricParameters.instrument.name=='VIRCAM') & \
(AstrometricParameters.filter.name=='ESO-Ks-0002') & \
(AstrometricParameters.quality_flags == 0) & \
(AstrometricParameters.is_valid >= 1) & \
(AstrometricParameters.creation_date > datetime.datetime(2017,5,14))).user_only()
sample=[]
for i in q :
    sample.append(i.reduced.filename)

lpu.run('Regrid', i='VIRCAM', red_filenames=sample, gra=cra, gdec=cdec, C=0)

# or
cra,cdec = 54.5 , -35.5
q = ((ReducedScienceFrame.instrument.name=='VIRCAM') & \
(ReducedScienceFrame.filter.name=='ESO-Ks-0002') &\
(ReducedScienceFrame.chip.name=='ESO-Virgo22') &\
(ReducedScienceFrame.creation_date>datetime.datetime(2017,11,1))).user_only()
N = len(q)
sample = list()
for i in range(0,N) :
    ap = (AstrometricParameters.reduced == q[i]).max('creation_date')
    if ap.RMS > 0:
        sample.append(q[i].filename)
        lpu.run('Regrid', i='VIRCAM', red_filenames=sample, gra=cra, gdec=cdec, C=1)

10. Coaddition

q2 = ((RegriddedFrame.is_valid==1) & \
(RegriddedFrame.creation_date > datetime.datetime(2017,11,23,15,31,17)) & \
(RegriddedFrame.quality_flags == 0) & \
(RegriddedFrame.is_valid >= 1) & \
(RegriddedFrame.filter.name=='ESO-Ks-0002') &\
(RegriddedFrame.chip.name=='ESO-Virgo22') &\
(RegriddedFrame.instrument.name=='VIRCAM'))
len(q2)
sample=[]
for i in q2 :
    sample.append(i.filename)
lpu.run('Coadd', instrument='VIRCAM',reg_filenames=sample , C=0 )

#############################################################
#############################################################

module unload awe
export PATH="/data/users/saifollahi/awe/Linux-redhat-7.3-x86_64/develop/astro/bin:$PATH"

#q = ReducedScienceFrame.select(instrument='VIRCAM', filter='ESO-Ks-0002')
q = ((ReducedScienceFrame.instrument.name=='VIRCAM') & \
(ReducedScienceFrame.filter.name=='ESO-Ks-0002') &\
(ReducedScienceFrame.creation_date>datetime.datetime(2018,3,1))).user_only()
N = len(q)
j=0
for i in range(0,N) :
    ap = (AstrometricParameters.reduced == q[i]).max('creation_date')
    if ap.RMS > 0:
        j=j+1
        #print (i, j, ap.RMS, ap.reduced.filename, ap.creation_date)
        #ap.reduced.retrieve()
        ap.reduced.weight.retrieve()
        reduced = str(ap.reduced.filename)
        #weight = str(ap.reduced.weight.filename)
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
        #command = 'mv ' + weight + ' ' + reduced + '.weight.fits'
        os.system(command)

#############################################################
#############################################################

q = ((ReducedScienceFrame.instrument.name=='VIRCAM') & \
(ReducedScienceFrame.filter.name=='ESO-Ks-0002') &\
(ReducedScienceFrame.creation_date>datetime.datetime(2017,11,1))).user_only()
N = len(q)
print (N)
j=0
for i in range(0,N) :
    ap = (AstrometricParameters.reduced == q[i]).max('creation_date')
    if ap.RMS < 0.001 and ap.quality_flags > 200 :
        #print (q[i].DATE_OBS)
        sample = list()
        sample.append(q[i].filename)
        #q[i].retrieve()
        #command='python header_corr.py '+ str(q[i].filename)
        #os.system(command)
        lpu.run('Astrometry', i='VIRCAM', red_filenames=sample, C=1)
        command='rm tmp*'
        os.system(command)


#############################################################
###################### D E L E T I N G ######################
#############################################################

q = ((AstrometricParameters.instrument.name=='VIRCAM') & \
(AstrometricParameters.filter.name=='ESO-Ks-0002') & \
(AstrometricParameters.creation_date > datetime.datetime(2017,11,8))).user_only()

for raw in q.project_only():
    context.delete(raw.inverse_objects(), delete_config=True, commit=True)
    context.delete(raw, delete_config=True, commit=True)

#############################################################

    q = ((ReducedScienceFrame.instrument.name == instrument_name) & \
    (ReducedScienceFrame.DATE_OBS < datetime.datetime(2012,12,1)) & \
    (ReducedScienceFrame.DATE_OBS > datetime.datetime(2012,1,1)) & \
    (ReducedScienceFrame.imstat.median<50) & \
    (ReducedScienceFrame.imstat.median>-10) & \
    (ReducedScienceFrame.creation_date > datetime.datetime(2018,3,1)) & \
    (ReducedScienceFrame.filter.name == filter_name)).user_only()
    N = len(q)
    print ('Number of total raw frames : ' + str(N))
    import os
    os.system('> swarp.list')
    for i in q :
        science = str(i.filename)
        print (science)
        i.retrieve()
        i.weight.retrieve()
        weight = str(i.weight.filename)
        os.system('mv '+weight+' '+science+'.weight.fits')
        print (science)
