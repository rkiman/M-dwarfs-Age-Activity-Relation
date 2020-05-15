#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.io.idl import readsav
#import matplotlib.pyplot as plt
import os
from astropy.table import Table

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

def calc_activity_fraction(color,ewha,method_bin):
    #polynomial for inactive stars
    #p_inactive = np.array([395.50786613,-3045.0560077,9964.56950987,
    #                       -18068.29454717,19816.23461458,
    #                       -13428.70983305,5469.77129071,
    #                       -1212.66808967,108.51204825])
    p_inactive = np.array([-22.81938807,108.69676223,-191.72732865,
                           148.87638458,-43.17475922])
    
    #define bins of colors 
    _,color_range = np.histogram(color,bins=method_bin)
    
    #Define a mask for active and inactive stars
    mask_inactive = abs(np.polyval(p_inactive,color)-ewha) < 0.75
    mask_active = ewha - np.polyval(p_inactive,color) >= 0.75
    
    n_bin = len(color_range)-1
    
    #arrays to be filled out as I calculate the active fraction
    color_bin = np.ones(n_bin)*np.nan
    active_fraction = np.ones(n_bin)*np.nan
    inactive_fraction = np.ones(n_bin)*np.nan
    n_array = np.ones(n_bin)*np.nan

    #calculation of the active fraction
    for i in range(n_bin):
        #select color range
        mask_color = (color_range[i] < color) * (color < color_range[i+1])
        n_tot = len(ewha[mask_color])
        if(n_tot!=0):
            active_fraction[i] = len(ewha[mask_color*mask_active])/n_tot
            inactive_fraction[i] = len(ewha[mask_color*mask_inactive])/n_tot
            color_bin[i] = np.nanmean(color[mask_color])
            n_array[i] = n_tot

    #calculate error of active fraction based on a binomial distribution
    error_binomial = (active_fraction*(1-active_fraction))/n_array

    return color_bin,active_fraction,np.sqrt(error_binomial),n_array

'''
def calc_activity_fraction(color,ewha,abs_mag,color_range,group=[],cluster='',
                           find_binaries=False):
    #polynomial for inactive stars
    
    p_inactive = np.array([395.50786613,-3045.0560077,9964.56950987,
                           -18068.29454717,19816.23461458,
                           -13428.70983305,5469.77129071,
                           -1212.66808967,108.51204825])
    
    mask_inactive = abs(np.polyval(p_inactive,color)-ewha) < 0.75
    mask_active = ewha - np.polyval(p_inactive,color) >= 0.75
    
    n_bin = len(color_range)-1
    
    if(find_binaries == False):
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
                
    elif(find_binaries == True):
        mask_binaries = find_binaries(color,abs_mag,group,cluster)

        color_bin = np.ones(n_bin)*np.nan
        active_fraction = np.ones(n_bin)*np.nan
        inactive_fraction = np.ones(n_bin)*np.nan
        n_array = np.ones(n_bin)*np.nan

        for i in range(n_bin):
            mask_color = (color_range[i] < color) * (color < color_range[i+1])
            mask_tot = mask_color*(~mask_binaries)
            n_tot = len(ewha[mask_tot])
            if(n_tot!=0):
                active_fraction[i] = len(ewha[mask_tot*mask_active])/n_tot
                inactive_fraction[i] = len(ewha[mask_tot*mask_inactive])/n_tot
                color_bin[i] = np.nanmean(color[mask_color])
                n_array[i] = n_tot

    return color_bin,active_fraction,inactive_fraction, n_array
'''