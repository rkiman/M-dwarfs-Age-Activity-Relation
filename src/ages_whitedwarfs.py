#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('/Users/rociokiman/Documents/wdwarfdate')
import wdwarfdate 

def calc_ages_wdm_binaries(wd_table):
    teff = wd_table['TeffH']
    e_teff = wd_table['e_TeffH']
    logg = wd_table['loggH']
    e_logg = wd_table['e_loggH']
    
    model_ifmr = 'Cummings_2018_MIST'
    
    results = wdwarfdate.calc_wd_age(teff,e_teff,logg,e_logg,
                                     method='bayesian',
                                     model_wd='DA',feh='p0.00',vvcrit='0.0',
                                     model_ifmr = model_ifmr,
                                     high_perc = 84, low_perc = 16,
                                     datatype='yr',
                                     plot=True,
                                     path='Catalogs/results_wd_ages/')
    
    return results


#wd_table = Table.read('../Catalogs/old/wdm_binaries.fits')

#print(wdwarfdate.__file__)
#calc_ages_wdm_binaries(wd_table)