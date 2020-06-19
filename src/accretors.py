#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from .astro import calc_number_single_stars

def identify_accretors(ls_compatible):
    '''
    Defines a mask_not_acc which indenty which stars are not accreators.
    '''
    
    N = len(ls_compatible)
    
    #Define mask to distiguis accreating stars according to the ewha
    mask_not_acc = np.ones(N)*np.nan
    spt = ls_compatible['spt']
    ewha = ls_compatible['ewha_all']
    mask_not_acc = []
    for x,y in zip(spt,ewha):
        mask_not_acc.append(def_mask_not_acc(x,y))
    mask_not_acc = np.array(mask_not_acc)
        
    #If one of the measurements of the same star is accreating, all are. 
    for i in range(1,int(max(ls_compatible['same_star']))+1):
        mask_same_star = ls_compatible['same_star'] == i
        if(any(mask_not_acc[mask_same_star]==False)):
            mask_not_acc[mask_same_star]=False
    
    log_file = open('log.txt','a')
    n_acc = len(ls_compatible[~mask_not_acc])
    log_file.write('Number of possible accretors: {}\n'.format(n_acc))
    n_single = calc_number_single_stars(ls_compatible[~mask_not_acc])
    text = 'Number of single possible accretors: {}\n'
    log_file.write(text.format(n_single))
    log_file.flush()
        
    return mask_not_acc


def spt_to_g_rp(spt):
    """
    Relation from Kiman et al. 2019
    """
    return -0.0036*spt**2 + 0.11*spt + 0.89

def def_mask_not_acc(spt,ewha):
    """
    White, R. J. & Basri, G. 
    VERY LOW MASS STARS AND BROWN DWARFS IN TAURUS-AURIGA. 
    Astrophys. J. 582, 1109–1122 (2003).
    """
    #If any is nan I cannot make a decision if they are accreating
    if(np.isnan(spt+ewha)):  
        return True
    else:
        if(spt >= -1 and spt<2.7):
            return ewha < 10.
        elif(spt >= 2.7 and spt<5.7):
            return ewha < 20.
        elif(spt >= 5.7 and spt<7.7):#elif(color>=spt_to_g_rp(5.7)):
            return ewha < 40.
        elif(spt>=7.7): #If it is later than M7.7 I don't have
            return True #any criteria
        else:
            return True
        
    
def calc_delta_ha_for_accretors(spt,ewha):
    """
    White, R. J. & Basri, G. 
    VERY LOW MASS STARS AND BROWN DWARFS IN TAURUS-AURIGA. 
    Astrophys. J. 582, 1109–1122 (2003).
    """
    #If any is nan I cannot make a decision if they are accreating
    if(np.isnan(spt+ewha)):
        return np.nan
    else:
        if(spt >= -1 and spt<2.7):
            return ewha - 10.
        elif(spt >= 2.7 and spt<5.7):
            return ewha - 20.
        elif(spt >= 5.7 and spt<7.7):#elif(color>=spt_to_g_rp(5.7)):
            return ewha - 40.
        else:
            return np.nan
        