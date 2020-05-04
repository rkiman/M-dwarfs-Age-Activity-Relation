#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt

def fit_inactive_sequence():

    path='/Users/rociokiman/Dropbox (Hunter College)/Catalogs/Newton2017_gaia.fit'
    newton2017 = Table.read(path)
    catalog = Table.read('../Catalogs/age_calibrators_bayes.fits')
    mask_pra = np.array(catalog['group_num'])==18
    catalog = catalog[mask_pra]
    
    color_newton = newton2017['phot_g_mean_mag'] - newton2017['phot_rp_mean_mag']
    color_calibrators = catalog['phot_g_mean_mag'] - catalog['phot_rp_mean_mag']
    
    ewha_newton = newton2017['EWHa']*(-1)
    ewha_calibrators = catalog['ewha']
    
    
    color = np.concatenate((np.array(color_newton),np.array(color_calibrators)))
    ewha = np.concatenate((np.array(ewha_newton),np.array(ewha_calibrators))) 
    
    mask_nan = ~np.isnan(color+ewha)* (ewha < 0.3)
    color1,ewha1 = color[mask_nan],ewha[mask_nan]
    
    p = np.polyfit(color1,ewha1,8)
    mask = abs(np.polyval(p,color1)-ewha1) < 0.2
    for i in range(10):
        p = np.polyfit(color1[mask],ewha1[mask],8)
        mask = abs(np.polyval(p,color1)-ewha1) < 0.2  
    
    x = np.linspace(0.5,1.5,10)
    
    plt.plot(color,ewha,'.')
    plt.plot(x,np.polyval(p,x),'-k')
    plt.ylim(-3,4)
    plt.xlim(0.4,1.6)
    plt.show()
    
    print(p)
    
fit_inactive_sequence()