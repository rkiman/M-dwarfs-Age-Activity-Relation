#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def check_gaia_match(compatibles):
    print('Checking Gaia cross-match')
    #Load columns I need
    same_star = compatibles['same_star']
    source_id = compatibles['source_id']
    #array to keep track of the index
    extra_id = np.arange(len(same_star))
    
    #Set parameters to save stars that have a bad crossmatch
    n_tot_diff_id = 0
    seen = []
    idx_to_remove = []
    
    for x in same_star:
        #Identify repeated stars
        mask_same_star = x == same_star
        if(len(same_star[mask_same_star])>1):
            #Make sure I don't look at the same star more than once
            if(x not in seen):
                seen.append(x)
                same_star_source_id = source_id[mask_same_star]
                #Check that all gaia ids are the same for all the stars
                if(any(same_star_source_id!=same_star_source_id[0])):
                    #If the gaia id is not the same for all the stars 
                    #identified as the same star, remove all of them
                    for idx_i in extra_id[mask_same_star]:
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
    same_star_good_match = compatibles_good_match['same_star']
    source_id_good_match = compatibles_good_match['source_id']
    n_tot_diff_id = 0
    seen = []
    for x in same_star_good_match:
        mask_same_star = x == same_star_good_match
        if(len(same_star_good_match[mask_same_star])>1):
            if(x not in seen):
                seen.append(x)
                same_star_source_id = source_id_good_match[mask_same_star]
                if(any(same_star_source_id!=same_star_source_id[0])):
                    n_tot_diff_id+=1
    assert n_tot_diff_id == 0
    
    return compatibles_good_match

