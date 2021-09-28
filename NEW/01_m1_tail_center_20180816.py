# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 09:44:53 2018

@author: IRLab
"""

binn = 1


import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits
from scipy.optimize import curve_fit
from glob import glob




def onclick(event):
    if event.button == 3: #right click
    #    global ix, iy
        ix, iy = event.xdata, event.ydata

        #test
        #ix, iy = 512, 543
                
 #       global sx, sy
        x = np.int(ix)
        y = np.int(iy)
        
        print 'x = %d, y = %d'%(x, y)

        
        # assign global variable to access outside of function
        global coords
        coords.append((x, y))
    
        # Disconnect after 2 clicks
        if len(coords) == 2:
            fig_full.canvas.mpl_disconnect(cid_full)
            plt.close(1)
            
        return

def func(x,a0,a1,a2,a3):
    return a0+a1*np.exp(-(x-a2)**2/(2*a3**2))

    
flist = glob('*.fits')
flist.sort()       
print '\nfile name = %s'%(flist[0])

## === get m0 coor
f = open('results/%s/00_m0_center_coor.txt'%(flist[0].split('.')[0]),'r')
line = f.readline()

cenY = float(line.split('\t')[0])
cenX = float(line.split('\t')[1])
FWHM = float(line.split('\t')[2])
hwidth = FWHM*1.5


print 'center y, x = %.2f, %.2f\n'%(cenY, cenX)

fitsucc = False
for u in range(0,3):
    if fitsucc == True:
        break
            
    coords = []    

    
    ## === show the full frame with center marked

    hdu = fits.open(flist[0])
    img, hdr = hdu[0].data, hdu[0].header
    
    z1, z2 = np.percentile(img,80), np.percentile(img,99)
    fig_full, ax_full = plt.subplots(figsize=(10,10))
    cax_full = ax_full.imshow(img, vmin=z1, vmax=z2, cmap='gray')
    fig_full.colorbar(cax_full)
    
    crosslen = 30
    ax_full.plot([cenX-crosslen/2, cenX+crosslen/2 ], [cenY, cenY], 'r-',lw = 1)
    ax_full.plot([cenX, cenX], [cenY-crosslen/2, cenY+crosslen/2], 'r-',lw = 1)
    
    cid_full = fig_full.canvas.mpl_connect('button_press_event', onclick)
    
    plt.show(1)
    
    ## ==== get line profile of m1 tail
    
    print coords
    sx, sy = coords[0]
        
    DLslope = (cenY-float(sy))/(cenX-float(sx))
    slope = -1.0/DLslope
    
    dx = hwidth/np.sqrt(1.0 + slope**2)
    dy = slope*dx
    
    x1, y1 = float(sx)+dx, float(sy)+dy
    x2, y2 = float(sx)-dx, float(sy)-dy
    
    
    print 'slope = %.2f\n'%(slope)
    
    profile = []
    
    if np.abs(slope) <= 1.0:
        if x1 <= x2:
            xi, yi = x1, y1
            xf = x2
        else:
            xi, yi = x2, y2
            xf = x1
        print 'xi, yi = %.2f, %.2f'%(xi, yi)        
        for x in range(int(xi), int(xf)):
            liney = round(slope*(x-sx) +sy)
            for y in range(int(liney-binn), int(liney+binn)):
                d = np.abs(slope*x-y -slope*sx+sy)/np.sqrt(slope**2+1)
                r = np.sqrt((x-xi)**2+(y-yi)**2 - d**2)
                if d <= binn/2.0:
  #                  print 'x, y, d, r, pix = %d, %d, %.2f, %.2f, %d \t O\n'%(x, y, d,r,img[y,x])
                    profile = profile + [(r, img[y,x])]
#                else:
 #                   print 'x, y, d, r, pix = %d, %d, %.2f, %.2f, %d \t X\n'%(x, y, d,r,img[y,x])
    if np.abs(slope) > 1.0:
        if y1 <= y2:
            xi, yi = x1, y1
            yf = y2
        else:
            xi, yi = x2, y2
            yf = y1
        print 'xi, yi = %.2f, %.2f'%(xi, yi)        
        for y in range(int(yi), int(yf)):
            linex = round((y-sy)/slope + sx)
            for x in range(int(linex-binn), int(linex+binn)):
                d = np.abs(slope*x-y -slope*sx+sy)/np.sqrt(slope**2+1)
                if d <= binn/2.0:
                    r = np.sqrt((x-xi)**2+(y-yi)**2-d**2)
#                    print 'x, y, d, r, pix = %d, %d, %.2f, %.2f, %d\n'%(x, y, d,r,img[y,x])
                    profile = profile + [(r, img[y,x])]
    
    profile.sort()
    profile = np.array(profile)
    print profile

    # === 1D gaussian fitting
    
    xdat, ydat = profile[:,0], profile[:,1]
    
    print ydat
    print np.max(ydat), np.min(ydat)
    print np.min(ydat), np.max(ydat)-np.min(ydat), np.median(xdat), FWHM/2.634
    try :
        popt, pcov = curve_fit(func, xdat, ydat, [np.min(ydat), np.max(ydat)-np.min(ydat), np.median(xdat), FWHM/2.634])
        print 'a =', popt[0]
        print 'b =', popt[1]
        print 'c =', popt[2]
        print 'd =', popt[3]
#        print pcov
    except RuntimeError :
        popt = [0.0, 0.0, 1.0, 0.0]
    
    xpl = np.linspace(0,int(hwidth*2),int(hwidth*2)*100)
    fig_plot2, ax_plot2 = plt.subplots(figsize=(7,7))
    ax_plot2.plot(xdat, ydat, 'bo',alpha=0.4)
    ax_plot2.plot(xpl, func(xpl, *popt), 'r-')
    fig_plot2.savefig('results/%s/01_1D_gaussian'%(flist[0].split('.')[0]))
    
    #plt.show(3)
    
    # === get m1 center
    
    if sx > xi:
        dx = (popt[2]-hwidth)/np.sqrt(1.0 + slope**2)
        dy = slope*dx
    else:
        dx = -(popt[2]-hwidth)/np.sqrt(1.0 + slope**2)
        dy = -slope*dx
    
    nx = sx + dx
    ny = sy + dy
    
    print 'x: %d->%.2f, y: %d->%.2f\n'%(sx,nx,sy,ny)
    
    
    DLslope = (cenY-ny)/(cenX-nx)
    #slope = -1.0/DLslope
    
    sx2, sy2 = coords[1]
    nx = nx + (sx2-sx)
    ny = ny + (sx2-sx)*DLslope    

    #full image with DL
    fig_DL, ax_DL = plt.subplots(figsize=(15,15))
    cax_DL = ax_DL.imshow(img, vmin=z1, vmax=z2, cmap='gray')
    ax_DL.plot([int(round(cenX)),int(round(nx))],[int(round(cenY)),int(round(ny))], 'r-', lw=1)
    ax_DL.plot([x1,x2],[y1,y2], 'g--', lw=1)
    fig_DL.colorbar(cax_DL)
    fig_DL.savefig('results/%s/01_dispersion_line'%(flist[0].split('.')[0]))
    
    plt.show(4)
    
    a = raw_input('stop?(y/n)')
    if a == 'y':
        fitsucc = True
    else:
        fitsucc = False
        
DLslope = (cenY-ny)/(cenX-nx)
angle = -np.arctan(DLslope)* 180.0 / np.pi
        
fout = open('results/%s/01_m1_center_coor.txt'%(flist[0].split('.')[0]), 'w')

if popt[0] != 0.0:
    fout.write('%.2f\t%.2f\t%.2f\n' % (ny, nx, angle))
    fitsucc = True        
else:
    fout.write('FITTING FAIL')
    fitsucc = False        

fout.close()   

