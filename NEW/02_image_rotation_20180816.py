# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 13:24:21 2018

@author: IRLab
"""

import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits
#from scipy.optimize import curve_fit
from glob import glob
#import os


flist = glob('*.fits')
flist.sort()       
print '\nfile name = %s'%(flist[0])
            
hdu = fits.open(flist[0])
img, hdr = hdu[0].data, hdu[0].header

img_sub = img.copy()#[400:650,400:650]

z1, z2 = np.percentile(img_sub,80), np.percentile(img_sub,99)
'''
fig_full, ax_full = plt.subplots(figsize=(10,10))
cax_full = ax_full.imshow(img_sub, vmin=z1, vmax=z2, cmap='gray')
fig_full.colorbar(cax_full)

plt.show(1)
'''

dtor = np.pi/180.0
lenth = len(img_sub)
print lenth




f = open('results/%s/01_m1_center_coor.txt'%(flist[0].split('.')[0]),'r')
line = f.readline()

m1y = float(line.split('\t')[0])
m1x =float(line.split('\t')[1])
theta = float(line.split('\t')[2])


f2 = open('results/%s/00_m0_center_coor.txt'%(flist[0].split('.')[0]),'r')
line2 = f2.readline()

m0y = float(line2.split('\t')[0])
m0x =float(line2.split('\t')[1])

#=========cal coor

fout = open('results/%s/02_coor_rotated.txt'%(flist[0].split('.')[0]), 'w')

m0x_rot = np.cos(theta*dtor)*m0x - np.sin(theta*dtor)*m0y + lenth*np.sin(theta*dtor)
m0y_rot = np.sin(theta*dtor)*m0x + np.cos(theta*dtor)*m0y

m1x_rot = np.cos(theta*dtor)*m1x - np.sin(theta*dtor)*m1y + lenth*np.sin(theta*dtor)
m1y_rot = np.sin(theta*dtor)*m1x + np.cos(theta*dtor)*m1y

fout.write('%.2f\t%.2f\n'% (m0y_rot, m0x_rot))
fout.write('%.2f\t%.2f\n'% (m1y_rot, m1x_rot))

   
fout.close()   

print('m0 coor = %.2f\t%.2f\n'% (m0y_rot, m0x_rot))
print('m1 coor = %.2f\t%.2f\n'% (m1y_rot, m1x_rot))

print lenth
lenth_rot = int(round(lenth*(np.cos(theta*dtor)+np.sin(theta*dtor))))
print lenth_rot
img_rot= np.zeros((lenth_rot+5, lenth_rot+5))
print img_rot

fxy = []
nxy = []

for y in range(0, lenth):
    print 'progress 1/2 %8d/%8d'%(y, lenth)
    for x in range(0, lenth):
        fx_rot = np.cos(theta*dtor)*x - np.sin(theta*dtor)*y + lenth*np.sin(theta*dtor)
        fy_rot = np.sin(theta*dtor)*x + np.cos(theta*dtor)*y
        nx_rot = int(round(fx_rot))
        ny_rot = int(round(fy_rot))
#        print 'x_rot, y_rot = %5d, %5d'%(nx_rot, ny_rot)
        #fxy = fxy + [(fx_rot, fy_rot)]
        #nxy = nxy + [(nx_rot, ny_rot)]
        if img_rot[ny_rot, nx_rot] != 0.0:
            img_rot[ny_rot, nx_rot] = (img_rot[ny_rot, nx_rot] + img_sub[y, x])/2.0
        else:
            img_rot[ny_rot, nx_rot] = img_sub[y, x]
        
for yr in range(0, lenth_rot):
    print 'progress 2/2 %8d/%8d'%(yr,lenth_rot)
    for xr in range(0, lenth_rot):
        if img_rot[yr, xr] == 0.0:
            img_rot[yr, xr] = (img_rot[yr-1, xr] + img_rot[yr, xr-1] + img_rot[yr + 1, xr] + img_rot[yr, xr + 1])/4.0

        
#fxy = np.array(fxy)
#nxy = np.array(nxy)

#fig, ax = plt.subplots(figsize=(10,10))
#ax.plot(fxy[:,0], fxy[:,1],'bo', alpha = 0.5)
#ax.plot(nxy[:,0], nxy[:,1],'g^', alpha = 0.5)

#plt.xticks(np.arange(0.5, 350.5, 1.0))
#plt.yticks(np.arange(0.5, 350.5, 1.0))
#plt.xlim(150, 160)
#plt.ylim(150, 160)

#plt.grid()


#plt.show(2)


fig, ax = plt.subplots(figsize=(10,10))
cax = ax.imshow(img_rot,vmin=z1, vmax=z2, cmap='gray')
fig.colorbar(cax)

fig1, ax1 = plt.subplots(figsize=(10,10))
cax1 = ax1.imshow(img_sub,vmin=z1, vmax=z2, cmap='gray')
fig1.colorbar(cax1)


plt.show(3)

fits.writeto('results/%s/test01.fits'%(flist[0].split('.')[0]), img_sub)
fits.writeto('results/%s/%s_rot.fits'%(flist[0].split('.')[0], flist[0].split('.')[0]), img_rot)

