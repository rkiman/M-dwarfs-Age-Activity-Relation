#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.io.idl import readsav
#import matplotlib.pyplot as plt
import os
from astropy.table import Table

def define_inactive_limit(color,ewha):
    mask_nan = ~np.isnan(color+ewha)* (ewha < 0.3)
    color,ewha = color[mask_nan],ewha[mask_nan]
    
    p = np.polyfit(color,ewha,3)
    mask = abs(np.polyval(p,color)-ewha) < 0.2
    for i in range(9):
        p = np.polyfit(color[mask],ewha[mask],8)
        mask = abs(np.polyval(p,color)-ewha) < 0.2  
    return p

def find_binaries(color,abs_mag,group,cluster):
    mask_cluster = group == cluster
    p=np.nan
    if(cluster == 'USCO'):
        return np.array([False for x in range(len(color))])
    elif(cluster == 'CARN'):
        return np.array([False for x in range(len(color))])
    elif(os.path.exists('Catalogs/Models/dr2_seq_fit_cf_'+cluster+'.sav')):
        model = readsav('Catalogs/Models/dr2_seq_fit_cf_'+cluster+'.sav')
        p = np.flipud(model['cf'])
    elif(cluster == 'PRA'):
        model = Table.read('Catalogs/Models/dr2_seq_fit_cf_PRA.fits')
        p = np.array(np.flipud(model['COEFFICIENT']))
    else:
        mask = ~np.isnan(color+abs_mag)*mask_cluster
        color1,abs_mag1 = color[mask],abs_mag[mask]
        p = np.polyfit(color1,abs_mag1,8)
        mask_outliers = abs((np.polyval(p,color1)- abs_mag1)) < 0.15
        for i in range(10):
            p = np.polyfit(color1[mask_outliers],abs_mag1[mask_outliers],8)
            mask_outliers = abs((np.polyval(p,color1)- abs_mag1)) < 0.15

    mask_binaries = abs((np.polyval(p,color)-0.75- abs_mag)) < 0.3
    '''
    x = np.linspace(np.nanmin(color[mask_cluster]),1.5,10)
    plt.scatter(color[mask_cluster],abs_mag[mask_cluster],s=10)
    plt.scatter(color[mask_binaries*mask_cluster],abs_mag[mask_binaries*mask_cluster],s=20,
                facecolors='none',edgecolor='k',linewidth=.5)
    plt.plot(x,np.polyval(p,x)-0.75,'--k',linewidth=0.5)
    plt.gca().invert_yaxis()
    plt.xlabel(r'$G-G_{\rm RP}$')
    plt.ylabel(r'$M_{\rm G}$')
    plt.savefig('results/cmd/'+cluster+'_CMD.png')
    plt.close()
    '''
    return mask_binaries

def calc_activity_fraction(color,ewha,abs_mag,group,cluster,color_range):
    #polynomial for inactive stars
    p = np.array([-373.08951374,2531.16938349,-7252.86298719,11428.40584469,
                  -10785.23563974,6209.01090379,-2121.21267524,398.27658421,
                  -34.62116487])
    
    n_bin = len(color_range)-1
    
    mask_inactive = abs(np.polyval(p,color)-ewha) < 1
    mask_active = (~mask_inactive) * (np.polyval(p,color) < ewha)
    
    if(np.logical_or(cluster == 'mlsdss',cluster == 'PRA')):
        color_bin = np.ones(n_bin)*np.nan
        active_fraction = np.ones(n_bin)*np.nan
        inactive_fraction = np.ones(n_bin)*np.nan
        n_array = np.ones(n_bin)*np.nan

        for i in range(n_bin):
            mask_color = (color_range[i] < color) * (color < color_range[i+1])
            mask_tot = mask_color
            n_tot = len(ewha[mask_tot])
            if(n_tot!=0):
                active_fraction[i] = len(ewha[mask_tot*mask_active])/n_tot
                inactive_fraction[i] = len(ewha[mask_tot*mask_inactive])/n_tot
                color_bin[i] = np.nanmean(color[mask_color])
                n_array[i] = n_tot
    else:
        mask_binaries = find_binaries(color,abs_mag,group,cluster)

        mask_cluster = group == cluster

        color_bin = np.ones(n_bin)*np.nan
        active_fraction = np.ones(n_bin)*np.nan
        inactive_fraction = np.ones(n_bin)*np.nan
        n_array = np.ones(n_bin)*np.nan

        for i in range(n_bin):
            mask_color = (color_range[i] < color) * (color < color_range[i+1])
            mask_tot = mask_color*mask_cluster*(~mask_binaries)
            n_tot = len(ewha[mask_tot])
            if(n_tot!=0):
                active_fraction[i] = len(ewha[mask_tot*mask_active])/n_tot
                inactive_fraction[i] = len(ewha[mask_tot*mask_inactive])/n_tot
                color_bin[i] = np.nanmean(color[mask_color])
                n_array[i] = n_tot

    return color_bin,active_fraction,inactive_fraction, n_array