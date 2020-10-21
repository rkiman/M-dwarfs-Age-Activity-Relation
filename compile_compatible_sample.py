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
star_index = src.find_repeated_stars(literature_search1['ra'],
                                    literature_search1['dec'])
star_index[star_index==0] = np.nan
literature_search1['star_index'] = star_index

#Remove haew nans, zeros and higher mass dwarfs from the sample. 
#We don't want them!
mask_nan = ~np.isnan(np.array(literature_search1['ewha']))
mask_zeros = np.array(literature_search1['ewha'])!=0
g = literature_search1['g_corr']
rp = literature_search1['rp_corr']
spt_column = src.get_spt(literature_search1['spt'], g-rp) #add photometric spt 
literature_search1['spt'] = spt_column                      
mask_m_dwarf = spt_column >= -1                           
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
#This also created the table summary of sources that is going to be used to 
#create the table for the paper
ls_compatible2 = src.select_compatible_measurements(literature_search,
                                                    max_order=2)
ls_compatible = src.check_gaia_match(ls_compatible2)

#Check number of repeated stars
seen = []
rep = 0
rep2 = 0
rep3 = 0
rep4 = 0
rep5 = 0
rep6 = 0
rep7 = 0
rep8 = 0
rep_more = 0
for x in ls_compatible['star_index'][~np.isnan(ls_compatible['star_index'])]:
    if(x not in seen):
        seen.append(x)
        mask_star_index = ls_compatible['star_index'] == x
        if(len(ls_compatible[mask_star_index])>1):
            rep += 1 
            if(len(ls_compatible[mask_star_index])==2):
                rep2+=1
            elif(len(ls_compatible[mask_star_index])==3):
                rep3+=1
            elif(len(ls_compatible[mask_star_index])==4):
                rep4+=1
            elif(len(ls_compatible[mask_star_index])==5):
                rep5+=1
            elif(len(ls_compatible[mask_star_index])==6):
                rep6+=1
            elif(len(ls_compatible[mask_star_index])==7):
                rep7+=1
            elif(len(ls_compatible[mask_star_index])==8):
                rep8+=1
            elif(len(ls_compatible[mask_star_index])>8):
                rep_more+=1
                
log_file.write('repeated_measurements:{}'.format(rep))
log_file.write('repeated_measurements 2:{}'.format(rep2))
log_file.write('repeated_measurements 3:{}'.format(rep3))
log_file.write('repeated_measurements 4:{}'.format(rep4))
log_file.write('repeated_measurements 5:{}'.format(rep5))
log_file.write('repeated_measurements 6:{}'.format(rep6))
log_file.write('repeated_measurements 7:{}'.format(rep7))
log_file.write('repeated_measurements 8:{}'.format(rep8))
log_file.write('repeated_measurements more:{}'.format(rep_more))

#Save sample
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
   
print('Done with compatible sample')
