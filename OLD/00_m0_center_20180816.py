# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

BOX = 20*2


import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits
from scipy.optimize import curve_fit
from glob import glob
import os




def func((x, y), a0, a1, a2, a3, a4):
    z= a0 + a1*np.exp(-((x - a2)**2+(y - a3)**2)/(2*a4**2))
    return z.ravel()


def onclick(event):
    if event.button == 3: #right click
        global ix, iy
        ix, iy = event.xdata, event.ydata

#        print 'x = %d, y = %d'%(ix, iy) #numpy.float64
        
        fig_full.canvas.mpl_disconnect(cid_full)
        plt.close(1)
        
        return



## === getting m0 coor by click
flist = glob('*.fits')
flist.sort()       
print '\nfile name = %s'%(flist[0])

try:
    os.mkdir('results/%s'%(flist[0].split('.')[0]))
except WindowsError:
    print ''
        
hdu = fits.open(flist[0])
img, hdr = hdu[0].data, hdu[0].header
    
    
fitsucc = False
for u in range(0,3):
    if fitsucc == True:
        break
        

    
    z1, z2 = np.percentile(img,80), np.percentile(img,99)
    fig_full, ax_full = plt.subplots(figsize=(5,5))
    cax_full = ax_full.imshow(img, vmin=z1, vmax=z2, cmap='gray')
    fig_full.colorbar(cax_full)
    
    cid_full = fig_full.canvas.mpl_connect('button_press_event', onclick)
    
    plt.show(1)
    
    ## === getting m0 box data
    
    sx = np.int(ix)
    sy = np.int(iy)
    #print 'sx, sy = %d, %d\n'%(sx, sy)
    
    subimg = img[(sy-BOX/2):(sy+BOX/2),(sx-BOX/2):(sx+BOX/2)]
    #        ax1.plot([sx-BOX/2,sx-BOX/2,sx+BOX/2,sx+BOX/2,sx-BOX/2], \
    #            [sy-BOX/2,sy+BOX/2,sy+BOX/2,sy-BOX/2,sy-BOX/2], 'r-',lw=2)
    
        
    ## ===  2D gaussian fitting 
    
    xdat = np.arange(0, BOX, 1.0)
    ydat = np.arange(0, BOX, 1.0)
    xdat, ydat = np.meshgrid(xdat, ydat)
    
    zdat = subimg.flatten()
    
    try :
        popt, pcov = curve_fit(func, (xdat, ydat), zdat)
            
        print 'a0 =', popt[0]
        print 'a1 =', popt[1]
        print 'a2 =', popt[2]
        print 'a3 =', popt[3]
        print 'a4 =', popt[4]
        
    except RuntimeError :
        popt = [0.0, 0.0, 0.0, 0.0, 1.0]
        print 'FITTING FAIL'
    
            
    ntvpl = 50
    xpl = np.linspace(0, 19, ntvpl)
    ypl = np.linspace(0, 19, ntvpl)
    xpl, ypl = np.meshgrid(xpl, ypl)
    zfit = func((xpl, ypl), *popt)
    #print zfit.reshape(ntvpl, ntvpl)
            
    fig_m0, ax_m0 = plt.subplots(figsize=(5,5))
                    #ax.hold(True)
    z1, z2 = np.percentile(zdat.reshape(BOX,BOX),30), np.percentile(zdat.reshape(BOX,BOX),100)
    cax_m0 = ax_m0.imshow(zdat.reshape(BOX,BOX), vmin=z1, vmax=z2, cmap='jet')
    ax_m0.contour(xpl, ypl, zfit.reshape(ntvpl,ntvpl), 5, colors='w')
    fig_m0.colorbar(cax_m0) 
    fig_m0.savefig('results/%s/00_contour'%(flist[0].split('.')[0]))
                    
    ## === center by 2D gaussian fitting
    
    cenX = sx - BOX/2 + popt[2]
    cenY = sy - BOX/2 + popt[3]
    FWHM = np.abs(popt[4]*2.35482)
    print 'center coor. x, y = %d -> %.2f ,%d -> %.2f \n'%(sx, cenX, sy, cenY)
    
    fout = open('results/%s/00_m0_center_coor.txt'%(flist[0].split('.')[0]), 'w') 
    if popt[0] != 0:     #when fitting fail (error)
        fout.write('%.2f\t%.2f\t%.2f\n' % (cenY, cenX, FWHM))
        fitsucc = True
    else:
        fout.write('FITTING FAIL')
        fitsucc = False
    fout.close()
            
    ## === radial plot
    radial_plot_data = []
    for i in range(0, 20):
        for j in range(0,20):
            radius = np.sqrt((i-popt[3])**2+(j-popt[2])**2)
    #        print '%.2f\t%6d\n' % (radius, subimg[i,j])
            radial_plot_data = radial_plot_data + [(radius, subimg[i,j])]
    radial_plot_data.sort()
    radial_plot_data = np.array(radial_plot_data)
    
    radial_plot_fit = []
    for i in range(0, 200):
        z= popt[0] + popt[1]*np.exp(-(i/10.0)**2/(2*popt[4]**2))
        radial_plot_fit = radial_plot_fit +[(i/10.0, z)]
    radial_plot_fit.sort()
    radial_plot_fit = np.array(radial_plot_fit)    
    
    fig_plot, ax_plot = plt.subplots(figsize=(7,7))
    ax_plot.plot(radial_plot_data[:,0], radial_plot_data[:,1], 'bo', alpha=0.4)
    ax_plot.plot(radial_plot_fit[:,0], radial_plot_fit[:,1] , 'r-')
    ax_plot.set_xlabel('Radius')
    ax_plot.set_ylabel('Pixel Value [ADU]')
    fig_plot.savefig('results/%s/00_radial_plot'%(flist[0].split('.')[0]))
        
    plt.show()
