#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.table import Table
import numpy as np
from .astro import add_corrected_magnitudes,get_spt,calc_number_single_stars
from .astro import calc_lhalbol
from .literature_search_compatible import find_repeated_stars
from .literature_search_compatible import select_compatible_measurements
from .check_gaia_match_func import check_gaia_match
from datetime import datetime
from .clear_gaia_kinematics import mask_gaia_cuts

def compile_compatible_sample():
    '''
    This function creates 
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
    literature_search1 = add_corrected_magnitudes(literature_search1)
    
    #Find repeated stars
    print('Finding repeated stars')
    star_index = find_repeated_stars(literature_search1['ra'],
                                     literature_search1['dec'])
    star_index[star_index==0] = np.nan
    literature_search1['star_index'] = star_index
    
    #Record total number of stars in the sample that have haew
    n_ls = len(literature_search1)
    text = 'Total number stars in the literature search sample: {}\n'
    log_file.write(text.format(n_ls))
    n_single = calc_number_single_stars(literature_search1)
    text = 'Total number of single stars: {}\n'
    log_file.write(text.format(n_single))
    
    #Select stars with good kinematics from gaia
    print('Making Gaia cuts')
    literature_search1 = mask_gaia_cuts(literature_search1)
    
    #Remove haew nans, zeros and higher mass dwarfs from the sample. 
    #We don't want them!
    mask_nan = ~np.isnan(np.array(literature_search1['ewha']))
    mask_zeros = np.array(literature_search1['ewha'])!=0
    g = literature_search1['g_corr']
    rp = literature_search1['rp_corr']
    spt_column = get_spt(literature_search1['spt'], g-rp)#add phot spt 
    literature_search1['spt'] = spt_column                      
    mask_m_dwarf = np.logical_or(spt_column >= -1, spt_column < 10)                          
    literature_search = literature_search1[mask_nan*mask_zeros*mask_m_dwarf]
    
    #Record number of M dwarfs in the sample that have haew
    n_ls = len(literature_search)
    text = 'Number M dwarfs with halpha in the literature search sample: {}\n'
    log_file.write(text.format(n_ls))
    n_single = calc_number_single_stars(literature_search)
    text = 'Number single M dwarfs with halpha: {}\n'
    log_file.write(text.format(n_single))
    
    #Select compatible catalogs. ls_compatible is the new literature search 
    #which has all the stars, but contains one column for haew only for the
    #compatible catalogs.
    #This also created the table summary of sources that is going to be used to 
    #create the table for the paper
    ls_compatible2 = select_compatible_measurements(literature_search,
                                                    max_order=2)
    ls_compatible = check_gaia_match(ls_compatible2)
    
    #Calculate LHalphaLbol
    print('Calculating Lhalbol')
    lhalbol,lhalbol_error = calc_lhalbol(ls_compatible['ewha'],
                                         ls_compatible['ewha_error'],
                                         ls_compatible['spt'])
    
    ls_compatible['lhalbol'] = lhalbol
    ls_compatible['lhalbol_error'] = lhalbol_error
    
    #Remove weird objects
    ra_weird = [279.9099730555555,279.90995833,279.909973,208.74074035]
                #23.808016111111108,21.11529166666666]
    dec_weird = [16.387083055555554,16.38711111,16.387083,5.21089823]
                # -7.2142930555555544,-33.91905555555555]
    mask_weird = []
    for rai,deci in zip(ls_compatible['ra'],ls_compatible['dec']):
        if(any(np.isclose(rai,ra_weird)) and any(np.isclose(deci,dec_weird))):
            mask_weird.append(True)
        else:
            mask_weird.append(False)
    mask_weird = np.array(mask_weird)
    
    ls_compatible = ls_compatible[~mask_weird]
    
    #Check number of repeated stars
    seen = []
    rep,rep2,rep3,rep4,rep5,rep6,rep7,rep8,rep_more = 0,0,0,0,0,0,0,0,0
    
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
                    
    log_file.write('repeated_measurements:{}\n'.format(rep))
    log_file.write('repeated_measurements 2:{}\n'.format(rep2))
    log_file.write('repeated_measurements 3:{}\n'.format(rep3))
    log_file.write('repeated_measurements 4:{}\n'.format(rep4))
    log_file.write('repeated_measurements 5:{}\n'.format(rep5))
    log_file.write('repeated_measurements 6:{}\n'.format(rep6))
    log_file.write('repeated_measurements 7:{}\n'.format(rep7))
    log_file.write('repeated_measurements 8:{}\n'.format(rep8))
    log_file.write('repeated_measurements more:{}\n'.format(rep_more))
    
    #Save sample
    ls_compatible.write('Catalogs/literature_search_gaia_compatible.fits',
                        format='fits')
    #Record number of compatible stars
    ls_compatible1 = ls_compatible[~np.isnan(ls_compatible['ewha'])]
    n_comp = len(ls_compatible1)
    text = 'Number of stars in the compatible sample: {}\n'
    log_file.write(text.format(n_comp))
    n_single = calc_number_single_stars(ls_compatible1)
    text = 'Number of single stars in the compatible sample: {}\n'
    log_file.write(text.format(n_single))
    log_file.flush()
       
    print('Done with compatible sample')
    return ls_compatible
