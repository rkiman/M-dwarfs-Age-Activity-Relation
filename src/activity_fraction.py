#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

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