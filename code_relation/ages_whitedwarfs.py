#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/rociokiman/Documents')
import wdwarfdate 

def calc_ages_wdm_binaries(wd_table):
    teff = wd_table['TeffH']
    e_teff = wd_table['e_TeffH']
    logg = wd_table['loggH']
    e_logg = wd_table['e_loggH']
    
    N = len(teff)
    results = np.ones((N,15))*np.nan
    model_ifmr = 'Cummings_2018_MIST'
    
    for i in range(N):
        results_i = wdwarfdate.calc_bayesian_wd_age(teff[i],e_teff[i], logg[i],
                                                     e_logg[i], n_mc=1000,
                                                     model_wd='DA', 
                                                     feh='p0.00',
                                                     vvcrit='0.0', 
                                                     model_ifmr = model_ifmr, 
                                                     n = 100, high_perc = 84, 
                                                     low_perc = 16, 
                                                     plot = True, 
                                                     save_dist = True, 
                                                     datatype = 'Gyr',
                                                     path = 'results/')
        results[i,:] = results_i
        
    table_name = 'Catalogs/wdm_binaries.fits'


    wd_table['ms_age_median'] = results[:,0]
    wd_table['ms_age_err_low'] = results[:,1]
    wd_table['ms_age_err_high'] = results[:,2]
    wd_table['cooling_age_median'] = results[:,3]
    wd_table['cooling_age_err_low'] = results[:,4]
    wd_table['cooling_age_err_high'] = results[:,5]
    wd_table['total_age_median'] = results[:,6]
    wd_table['total_age_err_low'] = results[:,7]
    wd_table['total_age_err_high'] = results[:,8]
    wd_table['initial_mass_median'] = results[:,9]
    wd_table['initial_mass_err_low'] = results[:,10]
    wd_table['initial_mass_err_high'] = results[:,11]
    wd_table['final_mass_median'] = results[:,12]
    wd_table['final_mass_err_low'] = results[:,13]
    wd_table['final_mass_err_high'] = results[:,14]

    #Save results to plot in notebook
    wd_table.write(table_name, format='fits', overwrite=True)
    
    return results


#wd_table = Table.read('../Catalogs/wdm_binaries.fits')

#calc_ages_wdm_binaries(wd_table)