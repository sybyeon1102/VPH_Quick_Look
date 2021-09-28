# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 15:17:39 2018

@author: Seoyeon Byeon
"""
datawidth = 10.0
binn = 4

import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits
from scipy.optimize import curve_fit
import os
from glob import glob


def func(x, a0, a1, a2, a3):
    return a0 + a1*np.exp(-(x-a2)**2/(2*a3**2))

get_pofiles = raw_input("get profiles?(y/n)\n")
get_fiterr = raw_input("get differences?(y/n)\n")

#print get_pofiles, get_fiterr

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

print flist[0].split('.')[0]

## ==== get line profile

#dispersion line slope    
DLslope = (cenY-tailY)/(cenX-tailX)
slope = -1.0/DLslope

inval= binn
deltX = float(inval)/np.sqrt(1.0 + DLslope**2)
if DLslope*(cenX-tailX) <0: 
    deltX = -deltX

hwidth = datawidth
dx = hwidth/np.sqrt(1.0 + slope**2)
dy = slope*dx

reNum = int(np.sqrt((cenX-tailX)**2+(cenY-tailY)**2)/inval +1.0)

fout = open('results/%s/02_par.txt'%(flist[0].split('.')[0]), 'w') 

chipl = []

for re in range(0, reNum):

    print 'progress...%d/%d'%(re+1,reNum)
    dlX, dlY = cenX + deltX*re , cenY + DLslope*deltX*re
    print 'dlX, dlY = %10.2f, %10.2f'%(dlX, dlY)
        
    x1, y1 = dlX + dx, dlY + dy
    x2, y2 = dlX - dx, dlY - dy
    
#    print 'slope = %.2f\n'%(slope)
    
    profile = []
    
    if np.abs(slope) <= 1.0:
        if x1 <= x2:
            xi, yi = x1, y1
            xf = x2
        else:
            xi, yi = x2, y2
            xf = x1
#        print 'xi, yi = %.2f, %.2f'%(xi, yi)        
        for x in range(int(xi), int(xf)):
            liney = round(slope*(x-dlX) + dlY)
            for y in range(int(liney-binn), int(liney+binn)):
                d = np.abs(slope*x-y -slope*dlX+dlY)/np.sqrt(slope**2+1)
                if d <= inval*0.5:                
                    r = np.sqrt((x-xi)**2+(y-yi)**2 - d**2)
 #                   print 'x, y, d, r, pix = %d, %d, %.2f, %.2f, %d \t O\n'%(x, y, d,r,img[y,x])
                    profile = profile + [(r, img[y,x])]

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
                    profile = profile + [(r, img[y,x])]
    
    profile.sort()
    profile = np.array(profile)

    plt.close()    
    
    
    # === 1D gaussian fitting
    
    xdat, ydat = profile[:,0], profile[:,1]
    
    try :
        popt, pcov = curve_fit(func, xdat, ydat)
#        print 'point number, bg, amp, center, sigma =%03d\t%10.2f\t%10.2f\t%10.2f\t%10.2f'%(np.size(xdat),popt[0], popt[1], popt[2], popt[3])

#        print pcov
    except RuntimeError :
        popt = [1.0, 0.0, 0.0, 1.0]
        
    ## for the unproper fitting
    #parameter
    '''
    
#0 fitting fail        popt[0] == 1.0
------------------------------------------
#1 small amp            popt[1]/popt[0] < 0.1 
------------------------------------------
#2 far center           abs(popt[2]-hwidth) > hwidth
#3 wide sigma           abs(popt[3]) < hwidth/3.0

0-1,2
1,2,3-1,2,3


    '''
    #parm
    if popt[0] == 1.0:
        parm = 1

    else:
        parm = 0
        if abs(popt[3]) > hwidth/3.0:
            parm = parm + 8
        if abs(popt[2]-hwidth) > 5.0:
            parm = parm + 4
        if popt[1]/popt[0] < 0.1:
            parm = parm + 2
    
    #center coords
    if np.abs(slope) <= 1.0:
        x_c = xi + popt[2] / np.sqrt(1.0 + slope**2)
        y_c = yi + popt[2] / np.sqrt(1.0 + slope**2) * slope
    else:
        y_c = xi + popt[2] / np.sqrt(1.0 + 1.0/slope**2)
        x_c = xi + popt[2] / np.sqrt(1.0 + 1.0/slope**2) / slope
    
        
    fout.write('%04d\t%10.2f\t%10.2f\t%10.2f\t%10.2f\t%d\t%10.2f\t%10.2f\n' % (re*binn, popt[0], popt[1], popt[2], popt[3], parm, x_c, y_c))
        


    #chi-sqr
    chiq = np.sum((ydat-func(xdat, *popt))**2)
    
    chipl = chipl + [(re*binn , np.log10(chiq), parm)]
#    print 're, chiq = %03d\t%20.2f\n'%(re, chiq)

    if get_pofiles == 'y' :#and re>= int((LL)/binn) and re <= int((UL )/binn):
        
        xpl = np.linspace(0,int(hwidth*2),int(hwidth*2)*100)
        fig_plot2, ax_plot2 = plt.subplots(figsize=(7,7))
        ax_plot2.plot(xdat, ydat, 'bo',alpha=0.4)
        ax_plot2.plot(xpl, func(xpl, *popt), 'r-')
        plt.ylim(np.min(ydat),np.max(ydat))
        plt.xlim(0,int(hwidth*2))
        
        try:
            os.mkdir('results/%s/02_1D_gaussian'%(flist[0].split('.')[0]))
        except WindowsError:
            print ''
        fig_plot2.savefig('results/%s/02_1D_gaussian/%04d'%(flist[0].split('.')[0],re*binn))

    if get_fiterr == 'y' :#and re>= LL and re <= UL:
        
        xpl = np.linspace(0,int(hwidth*2),int(hwidth*2)*100)
        fig_plot2, ax_plot2 = plt.subplots(figsize=(7,7))
        ax_plot2.plot(xdat, ydat-func(xdat, *popt), 'bo',alpha=0.4)
        ax_plot2.plot(xpl, np.zeros(np.size(xpl)), 'r-')
        z = np.max([-np.min(ydat-func(xdat, *popt)),np.max(ydat-func(xdat, *popt))])
        plt.ylim(-z,z)
        plt.xlim(0,int(hwidth*2))
        
        try:
            os.mkdir('results/%s/02_difference'%(flist[0].split('.')[0]))
        except WindowsError:
            print ''
        fig_plot2.savefig('results/%s/02_difference/%04d'%(flist[0].split('.')[0],re*binn))
    
    plt.close()
    

fout.close()

## ======chi2 graph
chipl = np.array(chipl)

fig_cpl, ax_cpl = plt.subplots(figsize=(20,10))
ax_cpl.plot(chipl[:,0], chipl[:,1],'bo', alpha = 0.4, label = 'data point')

plt.xlim(np.min(chipl[:,0])-1,np.max(chipl[:,0])+1)
plt.ylim(np.min(chipl[:,1])*0.9,np.max(chipl[:,1])*1.1 )
plt.grid()
'''
parmarr = chipl[:,2]


#fitting fail
index = np.where(parmarr != 1)
temp = chipl[:,1].copy()
temp[index] = 0

ax_cpl.plot(chipl[:,0], temp,'ko', label='fitting fail')

#wide sigma
index = np.where(parmarr < 8)
temp = chipl[:,1].copy()
temp[index] = 0

ax_cpl.plot(chipl[:,0], temp,'rs', label='wide sigma')

#far center
index = np.where(parmarr >= 8)
parmarr[index] = parmarr[index]-8

index = np.where(parmarr < 4)
temp = chipl[:,1].copy()
temp[index] = 0

ax_cpl.plot(chipl[:,0], temp,'c^', label='far center')

#low amp
index = np.where(parmarr >= 4)
parmarr[index] = parmarr[index]-4

index = np.where(parmarr < 2)
temp = chipl[:,1].copy()
temp[index] = 0

ax_cpl.plot(chipl[:,0],temp, 'w+', label='low amp')


plt.legend()
'''
plt.xticks(range(0,int(np.max(chipl[:,0])+1),10))


fig_cpl.savefig('results/%s/02_chi-sqr'%(flist[0].split('.')[0]))


plt.show()

fout = open('results/%s/02_index.txt'%(flist[0].split('.')[0]), 'w') 
getindex = True
while getindex == True:
    a = raw_input("input index (chi-sqr)\n")
    if a == '':
        getindex = False
    else:
        b= int(a)
        fout.write('%04d\n' % (b))

fout.close()

