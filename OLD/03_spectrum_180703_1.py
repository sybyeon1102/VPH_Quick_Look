# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 15:41:40 2018

@author: Seoyeon Byeon
"""
datawidth = 10.0
binn = 4

import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits
from glob import glob

flist = glob('*.fits')
flist.sort()       
print '\nfile name = %s'%(flist[0])

## === get m0 coor
f = open('results/%s/00_m0_center_coor.txt'%(flist[0].split('.')[0]),'r')
line = f.readline()

cenY= float(line.split('\t')[0])
cenX= float(line.split('\t')[1])

print 'center y, x = %.2f, %.2f\n'%(cenY, cenX)

## === get m1 tail coor
f = open('results/%s/01_m1_center_coor.txt'%(flist[0].split('.')[0]),'r')
line = f.readline()

tailY= float(line.split('\t')[0])
tailX= float(line.split('\t')[1])

print 'tail y, x = %.2f, %.2f\n'%(tailY, tailX)

## ==== get image array

hdu = fits.open(flist[0])

img, hdr = hdu[0].data, hdu[0].header



## ==== get line profile

#dispersion line slope    
DLslope = (cenY-tailY)/(cenX-tailX)
slope = -1.0/DLslope


inval= binn
deltX = inval/np.sqrt(1.0 + DLslope**2)
if DLslope >0:
    deltX = -deltX

reNum = int(np.sqrt((cenX-tailX)**2+(cenY-tailY)**2)/inval +1.0)

spectrum = []

dat = np.genfromtxt('results/%s/02_par.txt'%(flist[0].split('.')[0])) #re, BG, Amp, Center, STD, validity, xc, yc
dist = dat[:,0]
center = dat[:,3]
STD = dat[:,4]
validity = dat[:,5]
BGdat = dat[:,1]
BG = np.median(BGdat)
xc = dat[:,6]
yc = dat[:,7]

z1, z2 = np.percentile(img,80), np.percentile(img,99)
fig_DL, ax_DL = plt.subplots(figsize=(10,10))
cax_DL = ax_DL.imshow(img, vmin=z1, vmax=z2, cmap='gray')
fig_DL.colorbar(cax_DL)

fout = open('results/%s/03_ADU_pix.txt'%(flist[0].split('.')[0]), 'w')

for re in range(0, np.size(dist)):
    
#    dlX, dlY = cenX + deltX*re , cenY + DLslope*deltX*re
    dlX, dlY = xc[re], yc[re]
  #  print 'dlX, dlX2 = %10.2f\t%10.2f\t%10.2f\n'%(dlX, dlX2, dlX2-dlX)
 #   print '%10.2f\t%10.2f\n'%(dlX, dlY)
    
    hwidth = 3*np.abs(STD[re])

    if hwidth >datawidth:
        hwidth = datawidth
    if np.abs(center[re]-datawidth)>datawidth:
        dlX, dlY = cenX + deltX*re , cenY + DLslope*deltX*re
    '''
    hwidth = 10.0
    dlX, dlY = cenX + deltX*re , cenY + DLslope*deltX*re
    '''
    dx = hwidth/np.sqrt(1.0 + slope**2)
    dy = slope*dx
        
    x1, y1 = dlX + dx, dlY + dy
    x2, y2 = dlX - dx, dlY - dy
       
    totADU = 0.0
    
    if re/5.0-re/5 == 0 and validity[re] != 1:
        ax_DL.plot([x1,x2],[y1,y2], 'r-', lw=1)
    
    if np.abs(slope) <= 1.0:
        if x1 <= x2:
            xi, yi = x1, y1
            xf = x2
        else:
            xi, yi = x2, y2
            xf = x1
        for x in range(int(xi), int(xf)):
            liney = round(slope*(x-dlX) + dlY)
            for y in range(int(liney-binn), int(liney+binn)):
                d = np.abs(slope*x-y -slope*dlX+dlY)/np.sqrt(slope**2+1)
                r = np.sqrt((x-xi)**2+(y-yi)**2 - d**2)
                if d <= inval*0.5:
 #                   print 'x, y, d, r, pix = %d, %d, %.2f, %.2f, %d \t O\n'%(x, y, d,r,img[y,x])
                    totADU = totADU + img[y,x]-BG#[re]
  #              else:
  #                  print 'x, y, d, r, pix = %d, %d, %.2f, %.2f, %d \t X\n'%(x, y, d,r,img[y,x])
    if np.abs(slope) > 1.0:
        if y1 <= y2:
            xi, yi = x1, y1
            yf = y2
        else:
            xi, yi = x2, y2
            yf = y1
   #     print 'xi, yi = %.2f, %.2f'%(xi, yi)        
        for y in range(int(yi), int(yf)):
            linex = round((y-dlY)/slope + dlX)
            for x in range(int(linex-binn), int(linex+binn)):
                d = np.abs(slope*x-y -slope*dlX + dlY)/np.sqrt(slope**2+1)
                if d <= inval*0.5:
                    r = np.sqrt((x-xi)**2+(y-yi)**2-d**2)
   #                 print 'x, y, d, r, pix = %d, %d, %.2f, %.2f, %d\n'%(x, y, d,r,img[y,x])
                    totADU = totADU + img[y,x]-BG#[re]
 #   print BG#[re]
#    print totADU
    spectrum = spectrum + [(dist[re], totADU)]
    fout.write('%04d\t%10.2f\t%02d\n' % (dist[re], totADU/binn, validity[re]))
    
fout.close()

    
spectrum = np.array(spectrum)

C=1.0
D=0.0

#C=53.36
#D=-536.93

fig_sp, ax_sp = plt.subplots(figsize=(20,10))
ax_sp.plot(spectrum[:,0]*C+D, spectrum[:,1], 'bo-', alpha=0.4)

#fitting fail
index = np.where(validity != 1)
temp = spectrum[:,1].copy()
temp[index] = -10000

ax_sp.plot(spectrum[:,0]*C+D, temp,'ko', label='fitting fail')

#wide sigma
index = np.where(validity < 8)
temp = spectrum[:,1].copy()
temp[index] = -10000

ax_sp.plot(spectrum[:,0]*C+D, temp,'rs', label='wide sigma')

#far center
index = np.where(validity >= 8)
validity[index] = validity[index]-8

index = np.where(validity < 4)
temp = spectrum[:,1].copy()
temp[index] = -10000

ax_sp.plot(spectrum[:,0]*C+D, temp,'c^', label='far center')

#low amp
index = np.where(validity >= 4)
validity[index] = validity[index]-4

index = np.where(validity < 2)
temp = spectrum[:,1].copy()
temp[index] = -10000

ax_sp.plot(spectrum[:,0]*C+D,temp, 'w+', label='low amp')

#chi sqr
temp=np.zeros(np.size(spectrum[:,0]))-10000
temp2=spectrum[:,1].copy()

index = np.genfromtxt('results/%s/02_index.txt'%(flist[0].split('.')[0]))
if np.size(index)>1:
    for i in range(0, np.size(index)):
        re = int((index[i])/binn)
        temp[re] = temp2[re]
elif np.size(index) == 1:
    re = int((index)/binn)
    temp[re] = temp2[re]


ax_sp.plot(spectrum[:,0]*C+D, temp, 'kx', label='inaccurate fitting')

plt.legend()
plt.ylim(-5000,np.max(spectrum[:,1])*1.1)

plt.xticks(range(0,int(np.max(spectrum[:,0]*C+D)+1),10) )

fig_sp.savefig('results/%s/03_plot'%(flist[0].split('.')[0]))
fig_sp.savefig('results/03_plot_%s'%(flist[0].split('.')[0]))


fig_DL.savefig('results/%s/03_3sigma'%(flist[0].split('.')[0]))

plt.show()
