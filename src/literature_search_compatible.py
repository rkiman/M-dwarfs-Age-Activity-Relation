#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u

def select_compatible_measurements(literature_search,order):
    '''
    Select the catalogs which are first and second order compatible with
    Kiman et al. 2019. from the literature search sample.
    '''
    source_ref_table = Table.read('data/source_ref.csv')
    #Numbers for each reference
    source_num_ref = source_ref_table['source_num']
    
    ra = np.array(literature_search['ra'])
    dec = np.array(literature_search['dec'])
    ewha = np.array(literature_search['ewha'])
    ewha_err = np.array(literature_search['ewha_error'])
    source_num = np.array(literature_search['source_num'])
    
    #Identify repeated stars
    same_star = find_repeated_stars(ra,dec)
    
    #Identify repeated measurements
    mask_idx_to_remove = remove_repeated_measurements(ewha,ewha_err,same_star)
    
    #Using the repeated stars identify compatible catalogs
    mask_kiman = source_ref_table['source_ref'] == 'Kiman 2019'
    idx_kiman = source_ref_table['source_num'][mask_kiman][0]
    
    res = find_compatible_catalogs(idx_kiman,source_num_ref,same_star,
                                   source_num,ewha,order)
    compatible,prob,total_comp,order = res
    
    #Using the list of compatible catalogs make a mask of compatible stars
    #and define the final sample of compatible stars.
    mask_compatible = np.array([x in compatible for x in source_num])
    
    #Remove repeated measurements
    ls_compatible = literature_search[~mask_idx_to_remove]
    
    #Add a new ewha column to distinguish between compatible and not compatible
    ls_compatible['ewha_all'] = ls_compatible['ewha']
    ls_compatible['ewha_error_all'] = ls_compatible['ewha_error']
    ls_compatible['ewha'][~mask_compatible] = np.nan
    ls_compatible['ewha_error'][~mask_compatible] = np.nan
    
    ls_compatible.write('Catalogs/literature_search_gaia_compatible.fits',
                        format='fits')
    
    return ls_compatible
    
def find_repeated_stars(ra,dec):
    '''
    Finds repeated stars and assigns them the same number in the array 
    same_star. The number starts at 1. Stars with same_star = 0 are singles.
    The stars are consider the same if they are closer than 2arcsec.
    '''
    Ntot = len(ra)
    same_star = np.zeros(Ntot)
    c_all = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    dummy = 1
    print('Finding repeated stars:')
    for i in range(Ntot):
        if(i%10000==0):
            print('{} out of {}'.format(i,Ntot))
        c1 = SkyCoord(ra=ra[i]*u.deg, dec=dec[i]*u.deg)
        separation = c1.separation(c_all).arcsec
        mask_close = separation < 2
        if(len(ra[mask_close])>1 and any(same_star[mask_close]==0)):
            same_star[mask_close] = dummy
            dummy+=1
    return same_star

def remove_repeated_measurements(ewha,ewha_err,same_star):
    '''
    There are several stars repeated between the different catalogs, but
    some of the haew are exactly the same. This doesn't contribute to my 
    analyzis, so I remove them using this function.
    '''
    idx = np.arange(len(ewha))
    idx_to_remove = []
    for x in range(1,int(max(same_star)+1)):
        #Select group of the same star with the mask
        mask = (same_star == x) 
        if(len(ewha[mask])>1):
            #Select the ewha of the repeated star
            ewha_mask = np.round(ewha[mask],3)
            ewha_err_mask = ewha_err[mask]
            idx_mask = idx[mask]
            if(len(ewha_mask) != len(set(ewha_mask))):
                seen = []
                for x in ewha_mask:
                    mask_zero = ewha_mask-x == 0
                    if(len(ewha_mask[mask_zero])>1 and x not in seen):
                        idx_mask_zero = idx_mask[mask_zero]
                        mask_not_nan = ~np.isnan(ewha_err_mask[mask_zero])
                        if(len(idx_mask_zero[mask_not_nan])!=0):
                            idx_to_keep = idx_mask_zero[mask_not_nan][0]
                        else:
                            idx_to_keep = idx_mask_zero[0]
                        seen.append(x)
                        for y in idx_mask_zero[idx_mask_zero!=idx_to_keep]:
                            idx_to_remove.append(y)
                            
    mask_idx_to_remove = np.array([True if x in idx_to_remove 
                                   else False for x in idx])
    return mask_idx_to_remove

def find_compatible_catalogs(idx_kiman,source_num_ref,same_star,source_num,
                             ewha,order):
    
    res = calc_compatible_matrix(source_num_ref,same_star,source_num,ewha)
    matrix_prob,matrix_prob_all,ewha_i,ewha_j = res
    
    N_ref = len(source_num_ref)
    
    compatible = []
    prob = []
    total_comp = []
    order = []
    
    #First order compatible catalogs
    for i in range(N_ref):
        if(matrix_prob[i,idx_kiman] > 0.0):
            prob_i = matrix_prob[i,idx_kiman]/matrix_prob_all[i,idx_kiman]
            if(i not in compatible and prob_i > 0.9):
                compatible.append(i)
                prob.append(prob_i)
                total_comp.append(matrix_prob[i,idx_kiman])
                order.append(1)
    
    if(order==2):
        #Second order compatible catalogs
        for i in compatible:
            for j in range(N_ref):
                if((matrix_prob[j,i] > 0)):
                    prob_j = matrix_prob[j,i]/matrix_prob_all[j,i]
                    if(j not in compatible and prob_j > 0.9):
                        compatible.append(j)
                        prob.append(prob_j)
                        total_comp.append(matrix_prob[j,i])
                        order.append(2)
    compatible,prob = np.array(compatible),np.array(prob)
    total_comp,order = np.array(total_comp),np.array(order)
    
    #Sort comaptible catalogs
    idx = np.argsort(compatible)
    compatible = compatible[idx]
    prob = prob[idx]
    total_comp = total_comp[idx]
    order = order[idx]
    
    return compatible,prob,total_comp,order


def calc_compatible_matrix(source_num_ref,same_star,source_num,ewha):
    '''
    Calculates two matrices: 
        -matrix_prob_all is a (N,N) matrix where N is the number of catalogs 
        included in the literature search. The element matrix_prob[i,j] 
        (and the element matrix_prob[j,i]) contains how many stars the 
        catalog i has in common with the catalog j. 
        -matrix_prob is a (N,N) matrix where N is the number of catalogs 
        included in the literature search. The element matrix_prob[i,j] 
        (and the element matrix_prob[j,i]) contains how many stars the 
        catalog i has in common with the catalog j, which has a difference
        in halpha ew smaller than 3A. 
        
    '''
    N_ref = len(source_num_ref)
    ewha_i = []
    ewha_j = []

    matrix_prob = np.ones((N_ref,N_ref))*0
    matrix_prob_all = np.ones((N_ref,N_ref))*0
    
    for x in range(1,int(max(same_star)+1)):
        #Select group of the same star with the mask
        mask = (same_star == x) 
        if(len(ewha[mask])>1):            
            #Select the sources where the star is repeated
            source_num_mask = source_num[mask]
            #Select the ewha of the repeated star
            ewha_mask = ewha[mask]           
            N = len(source_num_mask)
            for i in range(N):
                for j in range(N):
                    idx_i = int(source_num_mask[i])
                    idx_j = int(source_num_mask[j])
                    if(idx_i!=idx_j):
                        #If the difference of the ewha is smaller than 3A,
                        #then add a 1 to matrix_prob
                        matrix_prob_all[idx_i,idx_j]+=1
                        if(abs(ewha_mask[i]-ewha_mask[j])<3):
                            matrix_prob[idx_i,idx_i]+=1
                    #Record the values of ewha for stars that are the same
                    #so they can be compared later. Add them only once (i<j).
                    if(idx_i<idx_j):
                        ewha_i.append(ewha_mask[i])
                        ewha_j.append(ewha_mask[j])
                        
    ewha_i, ewha_j = np.array(ewha_i), np.array(ewha_j)
    
    return matrix_prob,matrix_prob_all,ewha_i,ewha_j