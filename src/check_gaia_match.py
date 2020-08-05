#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def check_gaia_match(compatibles):
    print('Checking Gaia cross-match')
    #Load columns I need
    star_index = compatibles['star_index']
    source_id = compatibles['source_id']
    #array to keep track of the index
    extra_id = np.arange(len(star_index))
    
    #Set parameters to save stars that have a bad crossmatch
    n_tot_diff_id = 0
    seen = []
    idx_to_remove = []
    
    for x in star_index:
        #Identify repeated stars
        mask_star_index = x == star_index
        if(len(star_index[mask_star_index])>1):
            #Make sure I don't look at the same star more than once
            if(x not in seen):
                seen.append(x)
                star_index_source_id = source_id[mask_star_index]
                #Check that all gaia ids are the same for all the stars
                if(any(star_index_source_id!=star_index_source_id[0])):
                    #If the gaia id is not the same for all the stars 
                    #identified as the same star, remove all of them
                    for idx_i in extra_id[mask_star_index]:
                        idx_to_remove.append(idx_i)
                    n_tot_diff_id+=1
                    
    log_file = open('log.txt','a')
    text = 'Number of single stars removed because they had a bad cross-match\
    with gaia:{}\n'
    log_file.write(text.format(n_tot_diff_id))
    log_file.flush()
    #Make a mask to remove the stars with a bad cross-match
    mask_to_remove = np.array([x in idx_to_remove for x in extra_id])
    #New table with good cross-matchs
    compatibles_good_match = compatibles[~mask_to_remove]
    
    #Repeate code wit the new table to make sure that no one is left
    star_index_good_match = compatibles_good_match['star_index']
    source_id_good_match = compatibles_good_match['source_id']
    n_tot_diff_id = 0
    seen = []
    for x in star_index_good_match:
        mask_star_index = x == star_index_good_match
        if(len(star_index_good_match[mask_star_index])>1):
            if(x not in seen):
                seen.append(x)
                star_index_source_id = source_id_good_match[mask_star_index]
                if(any(star_index_source_id!=star_index_source_id[0])):
                    n_tot_diff_id+=1
    assert n_tot_diff_id == 0
    
    return compatibles_good_match

