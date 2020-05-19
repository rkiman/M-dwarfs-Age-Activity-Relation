#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.table import Table,vstack
import numpy as np
import os
import src
from datetime import datetime

'''
This main file runs code in the src folder to compile the sample of age
calibrators. First looks for m dwarfs in moving groups using Banyan Sigma. 
Then looks for white-dwarfs binaries and calculates ages for the white dwarfs.
The final sample is saved in a folder called Catalogs.

This main file creates (if they don't exists already, so delete before run if
they do):
    - 'Catalogs/literature_search_gaia_compatible.fits'
    - 'Catalogs/literature_search_accretors.fits'
    - 'Catalogs/literature_search_not_mg.fits'
    - 'Catalogs/literature_search_mg.fits'
    - 'Catalogs/wdm_binaries.fits'
    - 'Catalogs/age_calibrators_bayes.fits'
'''

#Open log file
log_file = open('log.txt','a')
log_file.write('\n')
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
log_file.write("Log's date and time: {}\n".format(dt_string))
log_file.flush()

###--------- Literature search sample ---------###

path_compatible = 'Catalogs/literature_search_gaia_compatible.fits'
if(os.path.exists(path_compatible)):
    ls_compatible = Table.read(path_compatible)
else:
    print('No file literature_search_gaia_compatible.fits')
    
    path = 'Catalogs/literature_search_gaia.fits'
    literature_search1 = Table.read(path)
    
    #Remove haew nans, zeros and higher mass dwarfs from the sample. 
    #We don't want them!
    mask_nan = ~np.isnan(np.array(literature_search1['ewha']))
    mask_zeros = np.array(literature_search1['ewha'])!=0
    g = literature_search1['phot_g_mean_mag']
    rp = literature_search1['phot_rp_mean_mag']
    mask_m_dwarf = g-rp > 0.8
    
    literature_search = literature_search1[mask_nan*mask_zeros*mask_m_dwarf]
    n_ls = len(literature_search)
    text = 'Number of stars in the literature search sample: {}\n'
    log_file.write(text.format(n_ls))
    
    #Add extinction corrected magnitudes to the sample
    literature_search = src.add_corrected_magnitudes(literature_search)
   
    #Find repeated stars
    same_star = src.find_repeated_stars(literature_search['ra'],
                                        literature_search['dec'])
    #same_star == 0 means that it doesn't have a repeated star.
    same_star[same_star==0] = np.nan
    literature_search['same_star'] = same_star
    
    #Select compatible catalogs
    ls_compatible = src.select_compatible_measurements(literature_search,
                                                       same_star,max_order=2)
    n_comp = len(ls_compatible[~np.isnan(ls_compatible['ewha'])])
    text = 'Number of stars in the compatible sample: {}\n'
    log_file.write(text.format(n_comp))
    log_file.flush()

#Identify possible accretors
mask_not_acc = src.identify_accretors(ls_compatible)
ls_c_not_acc = ls_compatible[mask_not_acc]
ls_c_acc = ls_compatible[~mask_not_acc]
ls_c_acc.write('Catalogs/literature_search_accretors.fits',format='fits',
               overwrite=True)

n_not_acc = len(ls_c_not_acc[~np.isnan(ls_c_not_acc['ewha'])])
log_file.write('Number of stars not accreating: {}\n'.format(n_not_acc))
log_file.flush()

#Identify M-dwarfs in moving groups
m_dwarfs_mg,m_dwarfs_not_mg = src.compile_m_moving_groups_sample(ls_c_not_acc)
log_file.flush()

#Identify M-dwarfs co-moving with a White dwarf
m_dwarfs_wd = src.compile_m_wd_sample(m_dwarfs_not_mg)
log_file.flush()

#Make Final Age-calibrators sample
#Combine results 
age_calibrators = vstack([m_dwarfs_mg,m_dwarfs_wd])

mask_ha = ~np.isnan(age_calibrators['ewha'])
n_tot_cal = len(age_calibrators[mask_ha])
log_file.write('Total number of age calibrators: {}'.format(n_tot_cal))

#Calculate LHalphaLbol
lhalbol,lhalbol_error = src.calc_lhalbol(age_calibrators['ewha'],
                                         age_calibrators['ewha_error'],
                                         age_calibrators['g_corr']-
                                         age_calibrators['rp_corr'])

age_calibrators['lhalbol'] = lhalbol
age_calibrators['lhalbol_error'] = lhalbol_error

#Save sample
age_calibrators.write('Catalogs/age_calibrators_bayes.fits', overwrite=True)


