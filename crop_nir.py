field1 = [57.076811, -34.450556]
field2 = [57.076811, -35.450556]
field3 = [57.076811, -36.450556]
field4 = [55.849239, -33.450556]
field5 = [55.849239, -34.450556]
field6 = [55.849239, -35.450556]
field7 = [55.849239, -36.450556]
field8 = [55.849239, -37.450556]
field9 = [54.621667, -33.450556]
field10 = [54.621667, -34.450556]
field11 = [54.621667, -35.450556]
field12 = [54.621667, -36.450556]
field13 = [54.621667, -37.450556]
field14 = [53.394095, -33.450556]
field15 = [53.394095, -34.450556]
field16 = [53.394095, -35.450556]
field17 = [53.394095, -36.450556]
field18 = [53.394095, -37.450556]
field19 = [52.166523, -34.450556]
field20 = [52.166523, -35.450556]
field21 = [52.166523, -36.450556]
field22 = [52.166523, -37.450556]
field23 = [50.938951, -34.450556]
field24 = [50.938951, -35.450556]
field25 = [50.938951, -36.450556]
field26 = [50.938951, -37.450556]
field27 = [49.711379, -36.450556]
field28 = [49.711379, -37.450556]
field29 = [0.,0.]
field30 = [57.076811, -33.450556]
field31 = [52.166523, -33.450556]
field32 = [50.938951, -33.450556]
field33 = [57.076811, -37.450556]
fields_ra = list()
fields_dec = list()
fields_ra.append(field1[0])
fields_dec.append(field1[1])
fields_ra.append(field2[0])
fields_dec.append(field2[1])
fields_ra.append(field3[0])
fields_dec.append(field3[1])
fields_ra.append(field4[0])
fields_dec.append(field4[1])
fields_ra.append(field5[0])
fields_dec.append(field5[1])
fields_ra.append(field6[0])
fields_dec.append(field6[1])
fields_ra.append(field7[0])
fields_dec.append(field7[1])
fields_ra.append(field8[0])
fields_dec.append(field8[1])
fields_ra.append(field9[0])
fields_dec.append(field9[1])
fields_ra.append(field10[0])
fields_dec.append(field10[1])
fields_ra.append(field11[0])
fields_dec.append(field11[1])
fields_ra.append(field12[0])
fields_dec.append(field12[1])
fields_ra.append(field13[0])
fields_dec.append(field13[1])
fields_ra.append(field14[0])
fields_dec.append(field14[1])
fields_ra.append(field15[0])
fields_dec.append(field15[1])
fields_ra.append(field16[0])
fields_dec.append(field16[1])
fields_ra.append(field17[0])
fields_dec.append(field17[1])
fields_ra.append(field18[0])
fields_dec.append(field18[1])
fields_ra.append(field19[0])
fields_dec.append(field19[1])
fields_ra.append(field20[0])
fields_dec.append(field20[1])
fields_ra.append(field21[0])
fields_dec.append(field21[1])
fields_ra.append(field22[0])
fields_dec.append(field22[1])
fields_ra.append(field23[0])
fields_dec.append(field23[1])
fields_ra.append(field24[0])
fields_dec.append(field24[1])
fields_ra.append(field25[0])
fields_dec.append(field25[1])
fields_ra.append(field26[0])
fields_dec.append(field26[1])
fields_ra.append(field27[0])
fields_dec.append(field27[1])
fields_ra.append(field28[0])
fields_dec.append(field28[1])
fields_ra.append(field29[0])
fields_dec.append(field29[1])
fields_ra.append(field30[0])
fields_dec.append(field30[1])
fields_ra.append(field31[0])
fields_dec.append(field31[1])
fields_ra.append(field32[0])
fields_dec.append(field32[1])
fields_ra.append(field33[0])
fields_dec.append(field33[1])
print fields_ra
print fields_dec

import os, sys
#for i in range(0,33) :
for i in range(0,2) :
    #os.system('rm k' + str(i+1) + '.*')
    command = 'swarp Ks-final.fits -c default.swarp -IMAGEOUT_NAME k' + str(i+1)+'.fits -WEIGHTOUT_NAME k' + str(i+1) + '.weight.fits -WEIGHT_TYPE MAP_WEIGHT' \
    + ' -WEIGHT_IMAGE Ks-final.weight.fits -CENTER '+str(fields_ra[i])+','+str(fields_dec[i])
    print command
    os.system(command)
