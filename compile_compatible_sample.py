#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.table import Table
import numpy as np

import src
from datetime import datetime

'''

This main file creates (if they don't exists already, so delete before run if
they do):
    - 'Catalogs/literature_search_gaia_compatible.fits'
'''

#Open log file and record time and date of the run
log_file = open('log.txt','a')
log_file.write('\n')
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
log_file.write("Log's date and time: {}\n".format(dt_string))
log_file.flush()

path = 'Catalogs/literature_search_gaia.fits'
literature_search1 = Table.read(path)

#Add extinction corrected magnitudes G and RP to the sample 
print('Adding corrected magnitudes')
literature_search1 = src.add_corrected_magnitudes(literature_search1)

#Find repeated stars
same_star = src.find_repeated_stars(literature_search1['ra'],
                                    literature_search1['dec'])
#same_star == 0 means that it doesn't have a repeated star, so I replaced
#it with a nan
same_star[same_star==0] = np.nan
literature_search1['same_star'] = same_star

#Remove haew nans, zeros and higher mass dwarfs from the sample. 
#We don't want them!
mask_nan = ~np.isnan(np.array(literature_search1['ewha']))
mask_zeros = np.array(literature_search1['ewha'])!=0
g = literature_search1['g_corr']
rp = literature_search1['rp_corr']
mask_m_dwarf = g-rp > 0.8
literature_search = literature_search1[mask_nan*mask_zeros*mask_m_dwarf]
#Record number of M dwarfs in the sample that have haew
n_ls = len(literature_search)
text = 'Number of stars in the literature search sample: {}\n'
log_file.write(text.format(n_ls))
n_single = src.calc_number_single_stars(literature_search)
text = 'Number of single stars in the literature search sample: {}\n'
log_file.write(text.format(n_single))

#Select compatible catalogs. ls_compatible is the new literature search 
#which has all the stars, but contains one column for haew only for the
#compatible catalogs.
ls_compatible2 = src.select_compatible_measurements(literature_search,
                                                    max_order=2)
ls_compatible = src.check_gaia_match(ls_compatible2)

ls_compatible.write('Catalogs/literature_search_gaia_compatible.fits',
                    format='fits')
#Record number of compatible stars
ls_compatible1 = ls_compatible[~np.isnan(ls_compatible['ewha'])]
n_comp = len(ls_compatible1)
text = 'Number of stars in the compatible sample: {}\n'
log_file.write(text.format(n_comp))
n_single = src.calc_number_single_stars(ls_compatible1)
text = 'Number of single stars in the compatible sample: {}\n'
log_file.write(text.format(n_single))
log_file.flush()
   

