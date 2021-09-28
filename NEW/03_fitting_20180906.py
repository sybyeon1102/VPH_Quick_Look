# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 12:56:49 2018

@author: Seoyeon Byeon
"""
import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits
from glob import glob
#import os
from scipy.optimize import curve_fit


def func(x, a0, a1, a2, a3):
    return a0 + a1*np.exp(-(x-a2)**2/(2*a3**2))


flist = glob('*.fits')
flist.sort()       
print '\nfile name = %s'%(flist[0])


##======= getting rotated coordinate
f = open('results/%s/02_coor_rotated.txt'%(flist[0].split('.')[0]), 'r')

'''
fout.write('%.2f\t%.2f\n'% (m0y_rot, m0x_rot))
fout.write('%.2f\t%.2f\n'% (m1y_rot, m1x_rot))
'''

line1 = f.readline()
line2 = f.readline()
print line1
print line2

m0y = float(line1.split('\t')[0])
m0x = float(line1.split('\t')[1])

m1y = float(line2.split('\t')[0])
m1x  =float(line2.split('\t')[1])

##====== getting m0 FWHM
f2 = open('results/%s/00_m0_center_coor.txt'%(flist[0].split('.')[0]), 'r')

line = f2.readline()

FWHM = float(line.split('\t')[2])


##========= getting rotated image 
hdu = fits.open('results/%s/%s_rot.fits'%(flist[0].split('.')[0], flist[0].split('.')[0]))
img, hdr = hdu[0].data, hdu[0].header

z1, z2 = np.percentile(img,80), np.percentile(img,99)
fig, ax = plt.subplots(figsize=(10,10))
cax = ax.imshow(img,vmin=z1, vmax=z2, cmap='gray')
fig.colorbar(cax)
ax.plot([m0x, m1x], [m0y, m1y], 'r-')

plt.show(1)




fout = open('results/%s/03_fitting_parm.txt'%(flist[0].split('.')[0]), 'w')
# when x = 620
cenY = (m0y + m1y) /2.0
#X = 460

if m0x < m1x:
    order = 1
else:
    order = -1
    
for X in range(int(round(m0x)), int(round(m1x)), order):
    
    print ''
    print ''
    print '######### X = ', X,'############' 
    print ''
    
    ##==== getting proper data width
    
    for count in range(0, 4):
        plt.close()
        print 'FWHM = ',FWHM
        xdat = np.arange(int(round(cenY-1.5*FWHM)), int(round(cenY+1.5*FWHM)), 1.0)
    #    print xdat, np.size(xdat)
        ydat = img[int(round(cenY-1.5*FWHM)):int(round(cenY+1.5*FWHM)), X].copy()
    #    print ydat, np.size(ydat)
        
        try :
            popt, pcov = curve_fit(func, xdat, ydat, [np.min(ydat), np.max(ydat)-np.min(ydat), cenY, FWHM/2.634])
            print 'bg, amp, center, sigma =%.2f\t%.2f\t%.2f\t%.2f'%(popt[0], popt[1], popt[2], popt[3])
            if popt[3] > 2.0 * FWHM/2.35482:
                popt[3] = 2.0 * FWHM/2.35482
            if popt[3] < 10.0/2.35482/3.0:
                popt[3] = 10.0/2.35482/3.0
                
        
        #        print pcov
        except Exception as ex:
            print 'Error name :', ex

    
        xpl = np.arange(int(round(cenY-1.5*FWHM)), int(round(cenY+1.5*FWHM)), 0.2)
        fig, ax = plt.subplots(figsize=(7,7))
        ax.plot(xdat, ydat, 'bo',alpha=0.4)
        ax.plot(xpl, func(xpl, *popt), 'r-')
            
        FWHM = 2.35482*np.abs(popt[3])
    print 'FWHM = ',FWHM 
  #  plt.show()
    plt.close()
        
    
    ##==== sigma clipping
    
    xdat = np.arange(int(round(cenY-1.5*FWHM)), int(round(cenY+1.5*FWHM)), 1.0)
#    print xdat, np.size(xdat)
    ydat = img[int(round(cenY-1.5*FWHM)):int(round(cenY+1.5*FWHM)), X].copy()
#    print ydat, np.size(ydat)
    
    xdatc = xdat.copy()
    ydatc = ydat.copy()
    
    
    
    for i in range(0,3):
        masking = np.ones(np.size(xdat))
        #print masking
        try :
            popt, pcov = curve_fit(func, xdat, ydat, [np.min(ydat), np.max(ydat)-np.min(ydat), cenY, FWHM/2.634])
            print 'bg, amp, center, sigma =%.2f\t%.2f\t%.2f\t%.2f'%(popt[0], popt[1], popt[2], popt[3])
            if popt[3] > 2.0 * FWHM/2.35482:
                popt[3] = 2.0 * FWHM/2.35482
            if popt[3] < 10.0/2.35482/3.0:
                popt[3] = 10.0/2.35482/3.0
        
        #        print pcov
        except Exception as ex:
            print 'Error name :', ex
            popt = [0.0, 0.0, 0.0, 1.0]

        
        #print func(xdat, *popt) - ydat 
    #    print np.size(xdat)
        diff = np.abs(func(xdat, *popt) - ydat)
    #    print 'difference',diff
        var = np.sum((diff - np.mean(diff))**2)/np.size(diff)
    #    print var
        sigma = np.sqrt(var)
        print 'sigma = ', sigma
        
        index =  np.where(diff > 3.0*sigma)
     #   print index
        masking[index] = 0.0
     #   print 'masking = ',masking
        xdat = xdat.copy()[np.where(masking == 1.0)]
        ydat = ydat.copy()[np.where(masking == 1.0)]
    
        xpl = np.arange(int(round(cenY-1.5*FWHM)), int(round(cenY+1.5*FWHM)), 0.2)
        fig, ax = plt.subplots(figsize=(7,7))
        ax.plot(xdat, ydat, 'bo',alpha=0.4)
        ax.plot(xpl, func(xpl, *popt), 'r-')
   #     plt.show()
        plt.close()
        
    fout.write('%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n'% (X, popt[0], popt[1], popt[2], 2.35482*np.abs(popt[3])))
fout.close()



