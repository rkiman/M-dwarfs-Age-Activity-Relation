#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 20:21:09 2020

@author: rociokiman
"""

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

def cross_match_to_white_dwarfs(m_dwarfs):
    '''
    Takes a sample of M dwarfs that don't belong to moving groups and 
    cross_matches it with a sample of white dwarfs from 
    Gentile Fusillo, N. Pietro et al. MNRAS, 482, 4570–4591 (2019).
    Return the cross-matched samples.
    '''
    #Load white dwarfs catalog
    w_dwarfs = Table.read('/Users/rociokiman/Documents/M-dwarfs-Age-Activity-Relation/Catalogs/wd_sources/gaiawd.fit')
    
    #Define the two catalogs in astropy SkyCoord format with ra and dec
    m_dwarfs_radec = SkyCoord(ra=np.array(m_dwarfs['ra'])*u.degree, 
                               dec=np.array(m_dwarfs['dec'])*u.degree)
    w_dwarfs_radec = SkyCoord(ra=np.array(w_dwarfs['RA_ICRS'])*u.degree, 
                               dec=np.array(w_dwarfs['DE_ICRS'])*u.degree)
    
    #Cross-match the two samples with a 10arcmin radius
    idx_md, idx_wd, d2d,_ = w_dwarfs_radec.search_around_sky(m_dwarfs_radec, 
                                                              10*u.arcmin)
    
    #Transform separation to arcmin to check in the test file
    separation = d2d.to(u.arcmin)/u.arcmin
    
    return m_dwarfs[idx_md], w_dwarfs[idx_wd], separation

def get_mask_binaries(limits,params_m,params_wd,pwd,f_pwd):
    '''
    Depending on the parameters of the white dwarfs and the m dwarfs, defines
    if they are co-movers and creates a mask to distinguish them.
    '''
    #Get the parameters for m and white dwarfs
    parallax_m, parallax_error_m, pmra_m, pmra_error_m, pmdec_m, pmdec_error_m = params_m
    parallax_wd, parallax_error_wd, pmra_wd, pmra_error_wd, pmdec_wd, pmdec_error_wd = params_wd
    
    #Get the limits in charge of defining if the stars are co-moving
    a1,a2,b,c = limits

    #Mask to identify high likelihood WD 
    #(See Gentile Fusillo, N. Pietro et al. MNRAS, 482, 4570–4591 (2019)):
    mask_all = pwd > 0.75 
    mask_limit = f_pwd == 1 
    mask_wd = np.logical_or(mask_all,mask_limit) 
    
    #Select stars with good parallaxes
    mask_parallax = (parallax_m/parallax_error_m > a1)*(parallax_wd/parallax_error_wd > a2)
    
    #Select co-movers that have close proper motion
    mask_pm = ((abs(pmra_m-pmra_wd) < b*(pmra_error_m+pmra_error_wd))
               *(abs(pmdec_m-pmdec_wd) < b*(pmdec_error_m+pmdec_error_wd)))
    #Select co-movers that have close parallax
    mask_parallax_diff = abs(parallax_m-parallax_wd) < c*(parallax_error_m+parallax_error_wd)

    return mask_parallax*mask_pm*mask_parallax_diff*mask_wd

def randomize(a, b, c, d, Ntot_wd):
    
    # Generate the permutation index array.
    permutation = np.random.permutation(a.shape[0])
    # Shuffle the arrays by giving the permutation in the square brackets.
    shuffled_a = a[permutation]
    shuffled_b = b[permutation]
    shuffled_c = c[permutation]
    shuffled_d = d[permutation]
    
    return shuffled_a[:Ntot_wd], shuffled_b[:Ntot_wd], shuffled_c[:Ntot_wd], shuffled_d[:Ntot_wd]

def calc_pca(params_m_0, params_wd_0, pwd, f_pwd, limits):
    '''
    Calculates the probability of chance alignment for the white dwarfs m-dwarfs
    pairs by randomizing the proper motions of the white dwarfs. Positions of
    both m and white dwarfs remain unchanged. 
    '''
    
    #Load white dwarfs catalog
    w_dwarfs = Table.read('/Users/rociokiman/Documents/M-dwarfs-Age-Activity-Relation/Catalogs/wd_sources/gaiawd.fit')
    
    Ntot_wd = len(pwd)
    N = 100 #Number of repetition for the probability
    prob = np.zeros(Ntot_wd)
     
    #Set parameters of the white dwarfs original pairs
    parallax_wd,parallax_error_wd,pmra_wd,pmra_error_wd,pmdec_wd,pmdec_error_wd = params_wd_0

    for i in range(N):
        #Randomize proper motions and errors using the full catalog of white dwarfs
        pmra_wd_random1,pmra_error_wd_random1,pmdec_wd_random1,pmdec_error_wd_random1 = randomize(np.array(w_dwarfs['pmRA']), 
                                                                                                  np.array(w_dwarfs['e_pmRA']), 
                                                                                                  np.array(w_dwarfs['pmDE']), 
                                                                                                  np.array(w_dwarfs['e_pmDE']),
                                                                                                  Ntot_wd)
        #Define new parameters list with the random proper motions
        params_wd = [parallax_wd,parallax_error_wd, 
                     pmra_wd_random1,pmra_error_wd_random1, 
                     pmdec_wd_random1,pmdec_error_wd_random1]
        
        #Decide if the pairs are co-movers or not with the random proper motions
        mask = get_mask_binaries(limits, params_m_0, params_wd, pwd, f_pwd)

        #If there are any pairs, the probability of chance alignment increases 
        #for that m-dwarf
        prob[mask] += 1

    return np.array(prob)/N #Return normalized probability
