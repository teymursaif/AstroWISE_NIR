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
from scipy.interpolate import griddata
from matplotlib.colors import LogNorm

command = 'export PATH="/data/users/saifollahi/awe/Linux-redhat-7.3-x86_64/develop/astro/bin:$PATH"'
os.system(command)

directory = ''
#list_name = 'fits3.list'
#chip_name = 'Virgo35'

list_name = sys.argv[1]
chip_name = sys.argv[2]

fits_list = open(directory + list_name, 'r')
zp_cat = open('zp.cat', 'w')
zp = np.array([])
dzp = np.array([])
x = np.array([])
y = np.array([])
for fits_file_name in fits_list :
    print ('Reading FITS : ' + str(fits_file_name) + '\n')
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

    command0 = 'rm newimage.fits'
    os.system(command0)
    fits.writeto('newimage.fits')
    command1 = 'rm ' + str(directory+fits_name)
    command2 = 'mv newimage.fits ' + str(directory+fits_name)
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
        if i <= 24 :
            sex_output_filtered.write(lines)
        if i > 24 :
            line = lines.split()
            flag = int(line[18])
            fwhm = float(line[20])
            ell = float(line[22])
            if flag < 1 and fwhm > 1.0 and fwhm < 10.0 and ell < 0.1:
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
                    d = math.sqrt( (RA-ra)**2 + (DEC-dec)**2 )
                    if d*3600.0 < 0.5 and abs(m-M) < 50 and abs(dm+dM)< 10:
                        #print (ra, RA, dec, DEC, m-M, fwhm, flag, ell)
                        print (X, Y, m, dm, M, dM, m-M)
                        lzp = list()
                        ldzp = list()
                        lx = list()
                        ly = list()
                        lzp.append(M-m)
                        ldzp.append(math.sqrt(dm**2+dM**2))
                        lx.append(X)
                        ly.append(Y)
                        zp = np.append(zp,lzp)
                        dzp = np.append(dzp,ldzp)
                        x = np.append(x, lx)
                        y = np.append(y, ly)
                        zp_cat.write( str(X) + ' ' + str(Y) + ' ' + str(M) + ' ' + \
                        str(dM) + ' ' + str(m) + ' ' + str(dm) + ' ' + str(m-M) + \
                        ' ' + str(dM+dm) + ' ' + str(fwhm) + ' ' + str(ra) + ' ' \
                        + str(dec) + '\n' )
                        #regions.write('fk5; circle ' + str(RA) + ' ' + str(DEC) \
                        #+ ' 0.002d # text={zp=' + str(m-M) + '}\n')
            check = 1
        ref.close()
    regions.close()
    #print zp
    #n, bins, patches = plt.hist(zp, 20, normed=1, facecolor='green', alpha=0.75)
    sex_output.close()
    crossmatched.close()
    print ('\n-----------------------------------------------------------------\n')


ZP = abs(np.median(zp))
ZP_sigma = np.std(zp)
zp_cat.write('ZP ')
zp_cat.write(str(ZP) + ' ' + str(ZP_sigma) + '\n')
fits_list.close()
zp_cat.close()

# Reading ZP and ZP_sigma from zp.cat
zps = np.array([])
fwhms = np.array([])
x = np.array([])
y = np.array([])

zp_cat = open('zp.cat', 'r')
for lines in zp_cat :
    line = lines.split()
    if line[0] == 'ZP' :
        ZP = abs(float(line[1]))
        ZP_sigma = float(line[2])
zp_cat.close()

fig, axs = plt.subplots(nrows=2, ncols=3, sharex=False, figsize=(24, 16))
ax = axs[0,0]
ax2 = axs[0,1]
ax3 = axs[1,0]
ax4 = axs[1,1]

ax5 = axs[0,2]
ax6 = axs[1,2]

zp_cat = open('zp.cat', 'r')
for lines in zp_cat :
    line = lines.split()
    if line[0] != 'ZP' :
        X = float(line[0])
        Y = float(line[1])
        M = float(line[2])
        dM = float(line[3])
        m = float(line[4])
        dm = float(line[5])
        zp = abs(float(line[6]))
        dzp = float(line[7])
        fwhm = float(line[8])
        ra = float(line[9])
        dec = float(line[10])
        lzp = list()
        lfwhm = list()
        lx = list()
        ly = list()
        lzp.append(zp-ZP)
        lfwhm.append(fwhm)
        lx.append(X)
        ly.append(Y)
        if  dzp <= 0.05 and m > 12.0 and m < 16.0 : #abs(zp-ZP) <= 0.2 :
            print m, M+abs(zp), zp-ZP, dzp, zp, fwhm
            ax.errorbar( m, zp , yerr=dzp, fmt='.', c='b')
            ax6.scatter(X,Y, marker='o', c='b', s=2.5**(abs(m-10)), alpha=0.5)
            ax3.scatter(X,Y, marker='o', c='k', alpha=0.5)
            ax4.scatter(X,Y, marker='o', c='k', alpha=0.5)
            #ax5.errorbar( ra+dec , zp , yerr=dzp, fmt='.', c='g')
            zps = np.append(zps, lzp)
            fwhms = np.append(fwhms, lfwhm)
            x = np.append(x, lx)
            y = np.append(y, ly)

ZP_sigma = np.std(zps)
zp_cat.close()

command = 'mv zp.cat ' + str(chip_name) + '_zp.cat'
os.system(command)

ax.set_title('Zero Point Deviation')
ax.set_xlabel('Measured Magnitude (2MASS)')
ax.set_ylabel('Deviation from Zero-Point')
ax.plot([0, 20], [ZP, ZP], 'r-', lw=4)
ax.plot([0, 20], [ZP-ZP_sigma, ZP-ZP_sigma], 'k--', lw=2)
ax.plot([0, 20], [ZP+ZP_sigma, ZP+ZP_sigma], 'k--', lw=2)
ax.text(12.1, ZP+0.35, 'Median Zero point = ' + str(ZP), fontsize=10)
ax.text(12.1, ZP+0.3, 'ZP Standard Deviation = ' + str(ZP_sigma), fontsize=10)
ax.grid(True)
ax.set_xlim([12,14])
ax.set_ylim([ZP-0.5,ZP+0.5])

ax6.set_title('Positions of Stars')
ax6.set_xlim([0,2048])
ax6.set_ylim([0,2048])
ax6.set_xlabel('X')
ax6.set_ylabel('Y')
ax6.grid(True)

ax2.set_title('Zero Point Histogram')
n, bins, patches = ax2.hist(zps, 25, normed=1, facecolor='blue', alpha=0.5)
ax2.plot([0, 0], [0, 20], 'r-', lw=4)
ax2.plot([0-ZP_sigma, 0-ZP_sigma], [0, n.max()*2], 'k--', lw=2)
ax2.plot([0+ZP_sigma, 0+ZP_sigma], [0, n.max()*2], 'k--', lw=2)
ax2.set_xlim([0-3*ZP_sigma,0+3*ZP_sigma])
ax2.set_ylim([0,1.1*n.max()])
ax2.set_xlabel('Measured ZP deviations')
ax2.set_ylabel('Number')
ax2.grid(True)
#(mu, sigma) = norm.fit(zps)
#yy = mlab.normpdf( bins, mu, sigma)
#l = ax2.plot(bins, yy, 'k--', linewidth=2)

fwhm_mean = np.median(fwhms)
fwhm_std = np.std(fwhms)
ax5.set_title('FWHM Histogram')
n2, bins2, patches2 = ax5.hist(fwhms, 25, normed=1, facecolor='red', alpha=0.5)
ax5.plot([0, 0], [0, 20], 'r-', lw=4)
ax5.plot([0-ZP_sigma, 0-ZP_sigma], [0, n.max()*2], 'k--', lw=2)
ax5.plot([0+ZP_sigma, 0+ZP_sigma], [0, n.max()*2], 'k--', lw=2)
ax5.set_xlim([0,6])
ax5.set_ylim([0,1.1*n2.max()])
ax5.set_xlabel('FWHM(pixel)')
ax5.set_ylabel('Number')
ax5.grid(True)
#(mu2, sigma2) = norm.fit(fwhm)
#yy2 = mlab.normpdf( bins2, mu2, sigma2)
#l2 = ax5.plot(bins2, yy2, 'k--', linewidth=2)

grid_x, grid_y = np.mgrid[0:2048:1024j, 0:2048:1024j]

grid_z0 = griddata((x,y), fwhms, (grid_x, grid_y), method='linear')
#grid_z0 = griddata((x,y), fwhms, (grid_x, grid_y), method='nearest')
#grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
#grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic')
cax3 = ax3.imshow(grid_z0.T, extent=(0,2048,0,2048), origin='lower')
fig.colorbar(cax3, cmap=cm.jet, ax=ax3)

grid_z1 = griddata((x,y), zps, (grid_x, grid_y), method='linear')
#grid_z0 = griddata((x,y), fwhms, (grid_x, grid_y), method='nearest')
#grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
#grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic')
cax4 = ax4.imshow(grid_z1.T, extent=(0,2048,0,2048), origin='lower')
fig.colorbar(cax4, cmap=cm.jet, ax=ax4)


#xi, yi = np.linspace(x.min(), x.max(), 128), np.linspace(y.min(), y.max(), 128)
#xi, yi = np.meshgrid(xi, yi)

#rbf = scipy.interpolate.Rbf(x, y, fwhms, function='linear')
#zi = rbf(xi, yi)
#cax3 = ax3.imshow(zi, vmin=0, vmax=10, origin='lower',
#           extent=[x.min(), x.max(), y.min(), y.max()])
#cbar = fig.colorbar(cax3, cmap=cm.jet, ax=ax3, ticks=[1.0, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0])

#rbf = scipy.interpolate.Rbf(x, y, zps, function='linear')
#zi = rbf(xi, yi)
#cax4 = ax4.imshow(zi, cmap=cm.jet, vmin=-3*ZP_sigma, vmax=3*ZP_sigma, origin='lower',
#           extent=[x.min(), x.max(), y.min(), y.max()])
#cbar = fig.colorbar(cax4, ax=ax4, ticks=[-0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2])


"""
from matplotlib.patches import Rectangle

size_mesh = 128
n_mesh = (2048/size_mesh)
N_mesh = (2048/size_mesh)**2
for i in range(0,N_mesh) :
    x =
    y =
    fwhm_map = np.array([])
    zp_map = np.array([])
    zp_cat = open('zp.cat', 'r')
    for lines in zp_cat :
        line = lines.split()
        if line[0] != 'ZP' :
            X = float(line[0])
            Y = float(line[1])
            M = float(line[2])
            dM = float(line[3])
            m = float(line[4])
            dm = float(line[5])
            zp = abs(float(line[6]))
            dzp = float(line[7])
            fwhm = float(line[8])
            lzp = list()
            lfwhm = list()
            lzp.append(zp)
            lfwhm.append(fwhm)
            nx = int(X/size_mesh)
            ny = int(Y/size_mesh)
            xy = int(ny * n_mesh + nx)
            xy = int(ny * n_mesh + nx)
            if xy==i :
                zp_map = np.append(zp_map, lzp)
                fwhm_map = np.append(fwhm_map, lfwhm)
        zp_median = np.median(zp_map)
        fwhm_median = np.median(fwhm_map)
        zp_std = np.std(zp_map)
        fwhm_std = np.std(fwhm_map)
        ax3.
        ax4.
        ax6.


    zp_cat.close()
    #print zp_map
    #print fwhm_map
    #print "----------------------"
"""

ax3.set_title('FWHM 2D Distribution (Pixel)')
ax3.set_xlim([0,2048])
ax3.set_ylim([0,2048])
ax3.set_xlabel('X')
ax3.set_ylabel('Y')
ax3.grid(True)


ax4.set_title('Zero Point 2D Distribution')
ax4.set_xlim([0,2048])
ax4.set_ylim([0,2048])
ax4.set_xlabel('X')
ax4.set_ylabel('Y')
ax4.grid(True)

#ax6.set_title('STDDEV of Zero Point 2D distribution')
#ax6.set_xlim([0,2048])
#ax6.set_ylim([0,2048])
#ax6.set_xlabel('RA')
#ax6.set_ylabel('DEC')
#ax6.grid(True)

fig.suptitle('Statistics of VIRCAM ' + str(chip_name) + '\n')
plt.savefig(chip_name+'-statistics.png', dpi=500)
#plt.show()
plt.close()

#ZP-Magnitude Diagram
#Histogram
#n, bins, patches = plt.hist(zp, 20, normed=1, facecolor='green', alpha=0.75)
#FWHM distribution
#zeropoint distribution plot (2D)
