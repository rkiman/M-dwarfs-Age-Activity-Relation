#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 16:34:55 2020

@author: rociokiman
"""

from astropy.table import Table
import numpy as np


def calc_lhalbol(ewha,ewha_error,g_rp):
    chi_douglas2014 = np.array([6.6453, 6.0334, 5.2658, 4.4872, 3.5926, 
                                2.4768, 1.7363, 1.2057, 0.6122, 0.3522])*1e-5
    g_rp_kiman2019 = np.array([0.93, 1.01, 1.09, 1.16, 1.23, 1.32, 1.41, 
                               1.47, 1.57, 1.63])
    p = np.polyfit(g_rp_kiman2019,chi_douglas2014,4)

    N = len(ewha)
    lhalbol = np.ones(N)*np.nan 
    lhalbol_err = np.ones(N)*np.nan
    for i in range(N):
        if((0.8 <= g_rp[i]) and (g_rp[i] <=1.65)):
            if(~np.isnan(ewha_error[i]+ewha[i])):
                dist_ewha = np.random.normal(ewha[i],ewha_error[i],2000)
                dist_lhalbol = dist_ewha*np.polyval(p,g_rp[i])
                lhalbol[i] = np.nanmedian(dist_lhalbol)
                lhalbol_err[i] = np.nanstd(dist_lhalbol)
            elif(np.isnan(ewha_error[i]) and ~np.isnan(ewha[i])):
                lhalbol[i] = ewha[i]*np.polyval(p,g_rp[i])
    
    return lhalbol,lhalbol_err
            
def organize_table_format(columns):
    
    labels = ['ra', 'dec', 'gaia_source_id', 'ra_gaia', 'dec_gaia', 'pmra',
              'pmra_error', 'pmdec', 'pmdec_error', 'parallax', 
              'parallax_error', 'phot_g_mean_flux', 'phot_g_mean_flux_error',
              'phot_g_mean_mag', 'phot_rp_mean_flux', 
              'phot_rp_mean_flux_error', 'phot_rp_mean_mag', 
              'phot_bp_mean_flux', 'phot_bp_mean_flux_error', 
              'phot_bp_mean_mag', 'ewha', 'ewha_error', 'ewha_all', 
              'ewha_error_all', 'lhalbol', 'lhalbol_error', 'age', 
              'age_error_low', 'age_error_high', 'group_num', 'group_name',
              'source_num', 'source_ref']
    N = len(labels)
    
    organized_table = Table()
    for i in range(N):
        organized_table[labels[i]] = np.array(columns[i])

    return organized_table
    
    
    
    