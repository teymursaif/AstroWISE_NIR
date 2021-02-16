
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
from matplotlib import cm

command = 'export PATH="/data/users/saifollahi/awe/Linux-redhat-7.3-x86_64/develop/astro/bin:$PATH"'
os.system(command)

directory = ''
#list_name = 'fits3.list'
#chip_name = 'Virgo35'

list_name = sys.argv[1]
chip_name = sys.argv[2]

###
fits_list = open(directory + list_name, 'r')
zp_cat = open('zp.cat', 'w')

for fits_file_name in fits_list :
    zp_cat.write('# ' +fits_file_name)
    zp = np.array([])
    dzp = np.array([])
    x = np.array([])
    y = np.array([])

    print ('Reading FITS : ' + str(fits_file_name))
    fits_file_name = fits_file_name.split()
    fits_file_name = fits_file_name[0]
    #weight_file_name = fits_file_name + '.weight.fits'
    head_file_name = fits_file_name + '.head'
    command = 'cp ' + str(fits_file_name) + ' image.fits'
    os.system(command)
    command = 'cp ' + str(head_file_name) + ' image.head'
    os.system(command)
    fits_file_name = 'image.fits'
    head_file_name = 'image.head'
    # Correcting the header

    header_corr = open(head_file_name)
    for lines in header_corr :
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

    #print (cd11,cd12,cd21,cd22)

    fits_name = fits_file_name
    fits = pyfits.open(directory+fits_name)
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

    #command0 = 'rm newimage.fits'
    #os.system(command0)
    fits.writeto('newimage2.fits')
    command1 = 'rm ' + str(directory+fits_name)
    command2 = 'mv newimage2.fits ' + str(directory+fits_name)
    os.system(command1)
    os.system(command2)
    fits.close()
    header_corr.close()

    #reading the header, ra, dec, etc
    fits_file = pyfits.open(directory + str(fits_file_name))
    header = fits_file[0].header
    RA = float(header['CRVAL1'])
    DEC = float(header['CRVAL2'])
    #print ra, dec

    #reading REFRENCE catalogue (2MASS) and limiting objects based on coordinates
    main_ref = open(directory + 'irsa_catalog_search_results.tbl', 'r')
    ref = open(directory+'ref.tbl', 'w')
    i = 0
    for lines in main_ref :
        i = i+1
        if i >= 72 :
            line = lines.split()
            ra = float(line[0])
            dec = float(line[1])
            frame_size = 0.4 * 2048.0 / 3600.0
            if abs(ra-RA) < frame_size and abs(dec-DEC) < frame_size :
                #print ra, dec, abs(ra-RA), frame_size
                ref.write(str(lines))
        elif i < 72 :
            ref.write(str(lines))
    main_ref.close()
    ref.close()

    #running sextractor and photometry
    #command1 = 'cp ' + directory + str(weight_file_name) + ' ' + directory + 'weight.fits'
    command2 = 'sex ' + str(fits_file_name)
    #os.system(command1)
    os.system(command2)

    #filter stars from noise in sex_output
    sex_output = open('sex_output.cat', 'r')
    sex_output_filtered = open('sex_output_filtered.cat', 'w')
    i = 0
    for lines in sex_output :
        i = i+1
        if i <= 27 :
            sex_output_filtered.write(lines)
        if i > 27 :
            line = lines.split()
            flag = int(line[18])
            fwhm = float(line[20])
            ell = float(line[22])
            if flag < 4 and fwhm > 1.0 and fwhm < 10.0 and ell < 0.1:
                sex_output_filtered.write(lines)
    sex_output.close()
    sex_output_filtered.close()
    #crossmatching sex_output.cat and ref.tbl
    sex_output = open('sex_output_filtered.cat', 'r')
    crossmatched = open('crossmatched.cat', 'w')
    regions = open('regions.reg', 'w')

    i = 0
    check = 0

    for lines1 in sex_output :
        i = i+1
        if i >= 25 :
            line1 = lines1.split()
            RA = float(line1[10])
            DEC = float(line1[11])
            X = float(line1[8])
            Y = float(line1[9])
            flag = int(line1[18])
            fwhm = float(line1[20])
            ell = float(line1[22])
            M = float(line1[6])
            dM = float(line1[7])
            #print RA,DEC
            j = 0
            ref = open('irsa_catalog_search_results.tbl', 'r')
            for lines2 in ref :
                j = j+1
                if j >= 72 :
                    line2 = lines2.split()
                    ra = float(line2[0])
                    dec = float(line2[1])
                    m = float(line2[14])
                    dm = (line2[15])
                    if dm == 'null' :
                        dm = 99.99
                    else :
                        dm = float(dm)
                    #print ra,dec
                    e = math.sqrt((dm)**2 + (dM)**2)
                    d = math.sqrt( (RA-ra)**2 + (DEC-dec)**2 )

                    if d*3600.0 < 0.5 and abs(m-M) < 50 and abs(dm+dM)< 10 \
                    and e < 0.1:
                        print X,Y,m,M, m-M
                        lzp = list()
                        lzp.append(m-M)
                        zp = np.append(zp,lzp)
                        zp_cat.write( 'D' + ' ' + str(X) + ' ' + str(Y) + ' ' + str(M) + ' ' + \
                        str(dM) + ' ' + str(m) + ' ' + str(dm) + ' ' + str(m-M) + \
                        ' ' + str(e) + ' ' + str(fwhm) + ' ' + str(ra) + ' '\
                        + str(dec) + '\n')
                        #regions.write('fk5; circle ' + str(RA) + ' ' + str(DEC) \
                        #+ ' 0.002d # text={zp=' + str(m-M) + '}\n')
            check = 1
        ref.close()
    regions.close()
    #print zp
    #n, bins, patches = plt.hist(zp, 20, normed=1, facecolor='green', alpha=0.75)
    sex_output.close()
    crossmatched.close()

    ZP = abs(np.median(zp))
    ZP_SIGMA = np.std(zp)
    zp_cat.write('## ZP ')
    zp_cat.write(str(ZP))
    zp_cat.write(' ZP_SIGMA ')
    zp_cat.write(str(ZP_SIGMA))
    zp_cat.write('\n###\n')

fits_list.close()
zp_cat.close()


command = 'mv zp.cat ' + str(chip_name) + '_zp.cat'
os.system(command)

###
ZP = list()
ZP_SIGMA = list()

zp_cat = open(str(chip_name) + '_zp.cat', 'r')
for lines in zp_cat :
    line = lines.split()
    flag = line[0]
    if flag == '##' :
        if line[2] != 'nan' and line[2] != 0.0:
            ZP.append(float(line[2]))
            ZP_SIGMA.append(float(line[4]))
        elif line[2] == 'nan' or line[2] == 0.0:
            ZP.append(99)
            ZP_SIGMA.append(99)
###############

zp_cat = open(str(chip_name) + '_zp.cat', 'r')
fig, axs = plt.subplots(nrows=4, ncols=1, sharex=False, figsize=(24, 18))
ax1 = axs[0]
ax2 = axs[1]
ax3 = axs[2]
ax4 = axs[3]
x = np.array([])
y = np.array([])
r = np.array([])
dzps = np.array([])
i = 0
for lines in zp_cat :
    line = lines.split()
    flag = line[0]
    if flag == '#' :
        i = i+1
        zp0 = ZP[i-1]
        zp_sigma0 = ZP_SIGMA[i-1]
    if flag == 'D' :
        X = float(line[1])
        Y = float(line[2])
        M = float(line[3])
        dM = float(line[4])
        m = float(line[5])
        dm = float(line[6])
        zp = abs(float(line[7]))
        e = float(line[8])
        fwhm = float(line[9])
        ra = float(line[10])
        dec = float(line[11])
        dzp = zp - zp0
        if abs(dzp) <= 0.2 and m > 13.0 and m < 16.0 and e < 0.1 :
            lx = list()
            ly = list()
            ldzp = list()
            lx.append(X)
            ly.append(Y)
            ldzp.append(dzp)
            x = np.append(x, lx)
            y = np.append(y, ly)
            dzps = np.append(dzps, ldzp)
            #and abs(ra-53.7) < 0.05 and abs(dec+35.5) < 0.1:
            R = math.sqrt( (X-1024)**2 + (Y-1024)**2 )
            lr = list()
            lr.append(R)
            r = np.append(r, lr)
            #print X, Y, m, dzp, e
            ax1.errorbar(X, dzp , yerr=e, fmt='.', c='b', alpha=.2)
            ax2.errorbar(Y, dzp , yerr=e, fmt='.', c='b', alpha=.2)
            ax3.errorbar(R, dzp, yerr=e, fmt='.', c='b', alpha=.2)
        elif abs(dzp) < 0.5 and m > 10.0 and m < 16.0 and e < 0.1 :
            ax1.errorbar(X, dzp , yerr=e, fmt='.', c='r', alpha=.2)
            ax2.errorbar(Y, dzp , yerr=e, fmt='.', c='r', alpha=.2)
            ax3.errorbar(R, dzp, yerr=e, fmt='.', c='r', alpha=.2)


import numpy as np
import pylab as plt

total_bins = 16

# Sample data

bins = np.linspace(x.min(),x.max(), total_bins)
delta = bins[1]-bins[0]
idx  = np.digitize(x,bins)
running_median = [np.median(dzps[idx==k]) for k in range(total_bins)]

ax1.scatter(x,dzps,color='r',alpha=.2,s=2)
ax1.plot(bins-delta/2,running_median,'r--',lw=4,alpha=.8)
ax1.axis('tight')

running_std = [dzps[idx==k].std() for k in range(total_bins)]
ax1.errorbar(bins-delta/2,running_median,
              running_std,fmt=None, lw=4, ecolor='r')

###########

bins = np.linspace(y.min(),y.max(), total_bins)
delta = bins[1]-bins[0]
idx  = np.digitize(y,bins)
running_median = [np.median(dzps[idx==k]) for k in range(total_bins)]

ax2.scatter(y,dzps,color='r',alpha=.2,s=2)
ax2.plot(bins-delta/2,running_median,'r--',lw=4,alpha=.8)
ax2.axis('tight')

running_std = [dzps[idx==k].std() for k in range(total_bins)]
ax2.errorbar(bins-delta/2,running_median,
              running_std,fmt=None, lw=4, ecolor='r')

#########

bins = np.linspace(r.min(),r.max(), total_bins)
delta = bins[1]-bins[0]
idx  = np.digitize(r,bins)
running_median = [np.median(dzps[idx==k]) for k in range(total_bins)]

ax3.scatter(r,dzps,color='r',alpha=.2,s=2)
ax3.plot(bins-delta/2,running_median,'r--',lw=4,alpha=.8)
ax3.axis('tight')

running_std = [dzps[idx==k].std() for k in range(total_bins)]
ax3.errorbar(bins-delta/2,running_median,
              running_std,fmt=None, lw=4, ecolor='r')


ax1.set_xlim([0,2048])
ax1.set_ylim([-0.2,0.2])
ax2.set_xlim([0,2048])
ax2.set_ylim([-0.2,0.2])
ax3.set_xlim([0,1400])
ax3.set_ylim([-0.2,0.2])

ax1.grid(True)
ax2.grid(True)
ax1.set_xlabel('X')
ax1.set_ylabel('Delta ZP')
ax2.set_xlabel('Y')
ax2.set_ylabel('Delta ZP')
ax3.grid(True)
ax3.set_xlabel('R')
ax3.set_ylabel('Delta ZP')
#plt.show()
zp_cat.close()

ax4.set_title('Zero Point Histogram')
n, bins, patches = ax4.hist(dzps, 16, normed=1, facecolor='blue', alpha=0.5)
ax4.set_xlabel('Delta Zero Points')
ax4.set_ylabel('Number')
ax4.grid(True)

fig.suptitle('Statistics of VIRCAM ' + str(chip_name) + '\n')
plt.savefig(str(chip_name)+'-zp-1night.png', dpi=300)
#plt.show()
#plt.close()
