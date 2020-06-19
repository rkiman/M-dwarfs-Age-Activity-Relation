#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.table import Table,vstack
import numpy as np
import src
from datetime import datetime
import os

'''
This main file runs code in the src folder to compile the sample of age
calibrators. First looks for m dwarfs in moving groups using Banyan Sigma. 
Then looks for white-dwarfs binaries and calculates ages for the white dwarfs.
The final sample is saved in a folder called Catalogs.

This main file creates (if they don't exists already, so delete before run if
they do):
    - 'Catalogs/literature_search_accretors.fits'
    - 'Catalogs/literature_search_all_groups.fits'
    - 'Catalogs/literature_search_not_mg.fits'
    - 'Catalogs/literature_search_mg.fits'
    - 'Catalogs/wdm_binaries.fits'
    - 'Catalogs/age_calibrators_bayes.fits'
'''

#Open log file and record time and date of the run
log_file = open('log.txt','a')
log_file.write('\n')
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
log_file.write("Log's date and time: {}\n".format(dt_string))
log_file.flush()

#Import literature search sample. If the selection of the compatible catalogs
#was already done, it just imports that file and skips the selection.
path_compatible = 'Catalogs/literature_search_gaia_compatible.fits'

ls_compatible = Table.read(path_compatible)

#Identify possible accretors
print('Identifying accretors')
mask_not_acc = src.identify_accretors(ls_compatible)
ls_c_not_acc = ls_compatible[mask_not_acc]
ls_c_acc = ls_compatible[~mask_not_acc]

ls_c_acc.write('Catalogs/literature_search_accretors.fits',format='fits',
               overwrite=True)

ls_c_not_acc1 = ls_c_not_acc[~np.isnan(ls_c_not_acc['ewha'])]
n_not_acc = len(ls_c_not_acc1)
log_file.write('Number of stars not accreating: {}\n'.format(n_not_acc))
n_single = src.calc_number_single_stars(ls_c_not_acc1)
text = 'Number of single stars in the compatible sample: {}\n'
log_file.write(text.format(n_single))
log_file.flush()

#Identify M-dwarfs in moving groups
print('Identifying moving group members')
m_dwarfs_mg,m_dwarfs_not_mg = src.compile_m_moving_groups_sample(ls_c_not_acc)
log_file.flush()

mask_not_mg = ~np.isnan(m_dwarfs_not_mg['ewha'])
n_singles = src.calc_number_single_stars(m_dwarfs_not_mg[mask_not_mg])
text = 'Number of single stars not in moving groups: {}\n'
log_file.write(text.format(n_single))
log_file.flush()

#Remove stars close to PRA area because is very crouded and not good to find
#white dwarfs binaries
ref_groups = Table.read('data/moving_groups_ref.csv')
pra_num = ref_groups['group_num'][ref_groups['name']=='PRA'][0]
m_dwarfs_not_mg = m_dwarfs_not_mg[m_dwarfs_not_mg['group_num']!=pra_num]

mask_not_mg = ~np.isnan(m_dwarfs_not_mg['ewha'])
n_singles = src.calc_number_single_stars(m_dwarfs_not_mg[mask_not_mg])
text = 'Number of single stars not in moving groups without PRA: {}\n'
log_file.write(text.format(n_single))
log_file.flush()

if(~os.path.exists('Catalogs/literature_search_not_mg.fits')):
    print('Saving literature_search_not_mg.fits')
    m_dwarfs_not_mg.write('Catalogs/literature_search_not_mg.fits',
                          format='fits')
        
#Identify M-dwarfs co-moving with a White dwarf
m_dwarfs_wd = src.compile_m_wd_sample(m_dwarfs_not_mg)
log_file.flush()

#Make Final Age-calibrators sample
#Combine results 
age_calibrators = vstack([m_dwarfs_mg,m_dwarfs_wd])

mask_ha = ~np.isnan(age_calibrators['ewha'])
n_tot_cal = len(age_calibrators[mask_ha])
log_file.write('Total number of age calibrators: {}\n'.format(n_tot_cal))
n_single = src.calc_number_single_stars(age_calibrators[mask_ha])
text = 'Total number of single age calibrators: {}\n'
log_file.write(text.format(n_single))
log_file.flush()

#Calculate LHalphaLbol
print('Calculating Lhalbol')
lhalbol,lhalbol_error = src.calc_lhalbol(age_calibrators['ewha'],
                                         age_calibrators['ewha_error'],
                                         age_calibrators['spt'])

age_calibrators['lhalbol'] = lhalbol
age_calibrators['lhalbol_error'] = lhalbol_error

#Remove weird objects
ra_weird = [279.9099730555555,279.90995833,279.909973,208.74074035,
            23.808016111111108]
dec_weird = [16.387083055555554,16.38711111,16.387083,5.21089823,
             -7.2142930555555544]
mask_weird = []
for ra_i,dec_i in zip(age_calibrators['ra'],age_calibrators['dec']):
    if(any(np.isclose(ra_i,ra_weird)) and any(np.isclose(dec_i,dec_weird))):
        mask_weird.append(True)
    else:
        mask_weird.append(False)
mask_weird = np.array(mask_weird)

age_calibrators = age_calibrators[~mask_weird]

#Save sample
print('Saving age calibrators sample')
age_calibrators.write('Catalogs/age_calibrators_bayes_not_clean.fits', 
                      overwrite=True)
print('Done')

