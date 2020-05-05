#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def identify_accretors(ls_compatible):
    '''
    Defines a mask_not_acc which indenty which stars are not accreators.
    '''
    
    N = len(ls_compatible)
    
    #Define mask to distiguis accreating stars according to the ewha
    mask_not_acc = np.ones(N)*np.nan
    g_rp = ls_compatible['phot_g_mean_mag']-ls_compatible['phot_rp_mean_mag']
    ewha = ls_compatible['ewha_all']
    mask_not_acc = []
    for x,y in zip(g_rp,ewha):
        mask_not_acc.append(def_mask_not_acc(x,y))
    mask_not_acc = np.array(mask_not_acc)
        
    #If one of the measurements of the same star is accreating, all are. 
    for i in range(1,int(max(ls_compatible['same_star']))+1):
        mask_same_star = ls_compatible['same_star'] == i
        if(any(mask_not_acc[mask_same_star]==False)):
            mask_not_acc[mask_same_star]=False
    
    n_acc = len(ls_compatible[~mask_not_acc])
    print('Number of possible accretors: {}'.format(n_acc))
            
    return mask_not_acc


def spt_to_g_rp(spt):
    """
    Relation from Kiman et al. 2019
    """
    return -0.0036*spt**2 + 0.11*spt + 0.89

def def_mask_not_acc(color,ewha):
    """
    White, R. J. & Basri, G. 
    VERY LOW MASS STARS AND BROWN DWARFS IN TAURUS-AURIGA. 
    Astrophys. J. 582, 1109–1122 (2003).
    """
    #If any is nan I cannot make a decision if they are accreating
    if(~np.isnan(color+ewha)):       
        if(color<spt_to_g_rp(2.7)):
            return ewha < 10
        elif(color<spt_to_g_rp(5.7)):
            return ewha < 20
        elif(color<spt_to_g_rp(7.7)):#elif(color>=spt_to_g_rp(5.7)):
            return ewha < 40
        elif(color>=spt_to_g_rp(7.7)): #If it is later than M7.7 I don't have
            return True                #any criteria
        else:
            return False
    else:
        return True
