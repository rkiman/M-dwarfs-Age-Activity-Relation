#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from astropy.table import Table
import sys
sys.path.append('/Users/rociokiman/Documents/wdwarfdate')
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
                                                     datatype = 'yr',
                                                     path = 'results/')
        results[i,:] = results_i

    
    return results


#wd_table = Table.read('../Catalogs/old/wdm_binaries.fits')

#print(wdwarfdate.__file__)
#calc_ages_wdm_binaries(wd_table)