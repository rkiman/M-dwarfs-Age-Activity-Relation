#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

def fit_inactive_sequence():

    deg = 3
    target_url = ('https://zenodo.org/record/2636692/files/'
                  'MLSDSS-GaiaDR2_extended.fits?download=1')
    mlsdss = fits.open(target_url)
    
    r_ext = mlsdss[1].data['EXTINCTION'][:,2]
    z_ext = mlsdss[1].data['EXTINCTION'][:,4]
    ext = r_ext-z_ext
    mask_ext = ext < 0.1
    subred = (mlsdss[1].data['photometric_sample_subred'] == 1) * mask_ext 
    
    ewha = mlsdss[1].data['EWHA'][subred]
    color = (mlsdss[1].data['phot_g_mean_mag'][subred] 
             - mlsdss[1].data['phot_rp_mean_mag'][subred])
    
    mask_nan = ~np.isnan(color+ewha)* (ewha < 1)
    color1,ewha1 = color[mask_nan],ewha[mask_nan]
    
    p = np.polyfit(color1,ewha1,deg)
    #mask = abs(np.polyval(p,color1)-ewha1) < 0.2
    #for i in range(2):
    #    p = np.polyfit(color1[mask],ewha1[mask],deg)
    #    mask = abs(np.polyval(p,color1)-ewha1) < 0.2  
    
    x = np.linspace(0.8,1.5,10)
    
    plt.plot(color,ewha,'.')
    plt.plot(x,np.polyval(p,x),'-k')
    plt.ylim(-3,4)
    plt.xlim(0.8,1.5)
    plt.show()
    
    print(p)
    
fit_inactive_sequence()