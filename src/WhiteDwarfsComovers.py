#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import os
from scipy.interpolate import interp1d
from .ages_whitedwarfs import calc_ages_wdm_binaries
from .astro import organize_table_format

def compile_m_wd_sample(m_dwarfs_not_mg,
                        file_name_binaries='wdm_binaries.fits'):
    #Cross-mach m-dwarf sample with white dwarf sample to find all the pairs in
    #a 10 arcmin radius
    print('Compiling sample of M dwarfs-white dwarfs binaries')
    m_dwarfs_pairs,w_dwarfs_pairs,_=cross_match_to_white_dwarfs(m_dwarfs_not_mg)
    
    #Make space removing tables that are not used anymore
    del m_dwarfs_not_mg
    
    #Define the parameters of the m-white dwarfs pairs
    params_m = [np.array(m_dwarfs_pairs['parallax']), 
                np.array(m_dwarfs_pairs['parallax_error']),
                np.array(m_dwarfs_pairs['pmra']), 
                np.array(m_dwarfs_pairs['pmra_error']),
                np.array(m_dwarfs_pairs['pmdec']), 
                np.array(m_dwarfs_pairs['pmdec_error'])]
    
    params_wd = [np.array(w_dwarfs_pairs['Plx']), 
                 np.array(w_dwarfs_pairs['e_Plx']),
                 np.array(w_dwarfs_pairs['pmRA']), 
                 np.array(w_dwarfs_pairs['e_pmRA']),
                 np.array(w_dwarfs_pairs['pmDE']), 
                 np.array(w_dwarfs_pairs['e_pmDE'])]
    
    #Define limits to decide what is a co-mover:
    limits = [8, #parallax snr cut for m dwarfs
              4, #parallax snr cut for white dwarfs
              3, #sigma difference for pmra and pmdec
              3] #sigma difference for parallax
    
    #Select co-movers with a mask from the pairs
    mask_comovers_all = get_mask_binaries(limits,params_m,params_wd,
                                          w_dwarfs_pairs['Pwd'],
                                          w_dwarfs_pairs['f_Pwd'])
    
    #Calculate probability of chance alignment for each pair
    prob_chance_align = calc_pca(params_m, params_wd,
                                 w_dwarfs_pairs['Pwd'], 
                                 w_dwarfs_pairs['f_Pwd'],
                                 limits)
    
    #Final mask for pairs includes selected co-movers which have a low 
    #probability of chance alignment
    mask_comovers = mask_comovers_all * (prob_chance_align <= 0.01)
    mask_comp = ~np.isnan(m_dwarfs_pairs['ewha'])
    
    log_file = open('log.txt','a')
    
    info = 'Number of m-dwarfs white dwarfs pairs: {}\n'
    N = len(np.array(m_dwarfs_pairs['parallax'])[mask_comovers_all*mask_comp])
    log_file.write(info.format(N))
    
    info = 'Number of m-dwarfs with prob of chance alignment > 0.01: {}\n'
    mask_info = mask_comovers_all * (prob_chance_align > 0.01) * mask_comp
    N = len(np.array(m_dwarfs_pairs['parallax'])[mask_info])
    log_file.write(info.format(N))
    
    info = 'Number of m-dwarfs white dwarfs pairs with prob ' \
           'of chance alignment <= 0.01: {}\n'
    N = len(np.array(m_dwarfs_pairs['parallax'])[mask_comovers*mask_comp])
    log_file.write(info.format(N))
    
    #Remove nans to calculate ages in the next step.
    mask_nan_teff_logg = ~np.isnan(w_dwarfs_pairs['TeffH']+
                                   w_dwarfs_pairs['loggH'])
    
    info='Number of m-dwarfs white dwarfs pairs I will calculate age for: {}\n'
    N = len(np.array(m_dwarfs_pairs['parallax'])[mask_comovers*mask_nan_teff_logg])
    log_file.write(info.format(N))
    log_file.flush()
    
    #Define the new sample of m-dwarfs that have a white dwarf co-moving
    m_co_movers = m_dwarfs_pairs[mask_comovers*mask_nan_teff_logg]
    w_co_movers = w_dwarfs_pairs[mask_comovers*mask_nan_teff_logg]

    #Calculate total age, cooling age, main sequence age, initial mass 
    #and final mass for the white dwarfs
    print('Calculating ages for white dwarfs')
    result_w_ages = calc_ages_wdm_binaries(w_co_movers)

    w_co_movers['ms_age_median'] = result_w_ages['ms_age_median']
    w_co_movers['ms_age_err_low'] = result_w_ages['ms_age_err_low']
    w_co_movers['ms_age_err_high'] = result_w_ages['ms_age_err_high']
    w_co_movers['cooling_age_median'] = result_w_ages['cooling_age_median']
    w_co_movers['cooling_age_err_low'] = result_w_ages['cooling_age_err_low']
    w_co_movers['cooling_age_err_high'] = result_w_ages['cooling_age_err_high']
    w_co_movers['total_age_median'] = result_w_ages['total_age_median']
    w_co_movers['total_age_err_low'] = result_w_ages['total_age_err_low']
    w_co_movers['total_age_err_high'] = result_w_ages['total_age_err_high']
    w_co_movers['initial_mass_median'] = result_w_ages['initial_mass_median']
    w_co_movers['initial_mass_err_low'] = result_w_ages['initial_mass_err_low']
    w_co_movers['initial_mass_err_high']=result_w_ages['initial_mass_err_high']
    w_co_movers['final_mass_median'] = result_w_ages['final_mass_median'] 
    w_co_movers['final_mass_err_low'] = result_w_ages['final_mass_err_low']
    w_co_movers['final_mass_err_high'] = result_w_ages['final_mass_err_high']
    
    m_w_co_movers = hstack([w_co_movers,m_co_movers])
    
    good_ages = idenfy_good_wd_ages(m_w_co_movers)

    m_w_co_movers['good_ages'] = good_ages
    
    m_w_co_movers.write('Catalogs/'+file_name_binaries, format = 'fits', 
                        overwrite = True)
        
    #Organize table format for future steps
    N_final = len(m_co_movers)
    
    columns = [m_co_movers['ra'], m_co_movers['dec'], m_co_movers['spt'],
               m_co_movers['source_id'],
               m_co_movers['ra_x'], m_co_movers['dec_x'], m_co_movers['pmra'], 
               m_co_movers['pmra_error'], m_co_movers['pmdec'],
               m_co_movers['pmdec_error'], m_co_movers['parallax'],
               m_co_movers['parallax_error'], m_co_movers['phot_g_mean_flux'],
               m_co_movers['phot_g_mean_flux_error'], 
               m_co_movers['phot_g_mean_mag'],m_co_movers['phot_rp_mean_flux'],
               m_co_movers['phot_rp_mean_flux_error'], 
               m_co_movers['phot_rp_mean_mag'],
               m_co_movers['phot_bp_mean_flux'],
               m_co_movers['phot_bp_mean_flux_error'], 
               m_co_movers['phot_bp_mean_mag'], m_co_movers['g_corr'],
               m_co_movers['rp_corr'], m_co_movers['ewha'], 
               m_co_movers['ewha_error'], m_co_movers['ewha_all'], 
               m_co_movers['ewha_error_all'], m_co_movers['lhalbol'],
               m_co_movers['lhalbol_error'], 
               w_co_movers['total_age_median'],
               w_co_movers['total_age_err_low'], 
               w_co_movers['total_age_err_high'],
               np.zeros(N_final), np.array(['WD' for i in range(N_final)]),
               m_co_movers['star_index'],
               m_co_movers['source_num'], m_co_movers['source_ref']]

    m_co_movers_organized = organize_table_format(columns)
    
    return m_co_movers_organized[good_ages==1]

def test_wd_table(table,table1):
    assert len(table) == len(table1)
    for x in range(10):
        i = int(np.random.uniform(0,len(table)))
        assert table['RA_ICRS'][i] == table1['RA_ICRS'][i]
        assert table['DE_ICRS'][i] == table1['DE_ICRS'][i]

def cross_match_to_white_dwarfs(m_dwarfs):
    '''
    Takes a sample of M dwarfs that don't belong to moving groups and 
    cross_matches it with a sample of white dwarfs from 
    Gentile Fusillo, N. Pietro et al. MNRAS, 482, 4570–4591 (2019).
    Return the cross-matched samples.
    '''

    if(os.path.exists('Catalogs/wd_sources/gaiawd.fit')):
        w_dwarfs = Table.read('Catalogs/wd_sources/gaiawd.fit')
        ra_wd = w_dwarfs['RA_ICRS']
        dec_wd = w_dwarfs['DE_ICRS']
    else:
        url="ftp://cdsarc.u-strasbg.fr/pub/cats/J/MNRAS/482/4570/gaia2wd.dat.gz" 
        w_dwarfs = Table.read(url, readme="data/ReadMe", format="ascii.cds")
        ra_wd = w_dwarfs['RAdeg']
        dec_wd = w_dwarfs['DEdeg']

    #Define the two catalogs in astropy SkyCoord format with ra and dec
    w_dwarfs_radec = SkyCoord(ra=np.array(ra_wd)*u.degree, 
                               dec=np.array(dec_wd)*u.degree)
    
    m_dwarfs_radec = SkyCoord(ra=np.array(m_dwarfs['ra'])*u.degree, 
                              dec=np.array(m_dwarfs['dec'])*u.degree)
    
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
    parallax_m,e_parallax_m,pmra_m,e_pmra_m,pmdec_m,e_pmdec_m = params_m
    parallax_wd,e_parallax_wd,pmra_wd,e_pmra_wd,pmdec_wd,e_pmdec_wd = params_wd
    
    #Get the limits in charge of defining if the stars are co-moving
    a1,a2,b,c = limits

    #Mask to identify high likelihood WD 
    #(See Gentile Fusillo, N. Pietro et al. MNRAS, 482, 4570–4591 (2019)):
    mask_all = pwd > 0.75 
    mask_limit = f_pwd == 1 
    mask_wd = np.logical_or(mask_all,mask_limit) 
    
    #Select stars with good parallaxes
    mask_parallax_m = parallax_m/e_parallax_m > a1
    mask_parallax_wd = parallax_wd/e_parallax_wd > a2
    mask_parallax = mask_parallax_m * mask_parallax_wd
    
    #Select co-movers that have close proper motion
    mask_pm = ((abs(pmra_m - pmra_wd) < b * (e_pmra_m + e_pmra_wd))
               * (abs(pmdec_m - pmdec_wd) < b * (e_pmdec_m + e_pmdec_wd)))
    #Select co-movers that have close parallax
    mask_dist = abs(parallax_m-parallax_wd) < c*(e_parallax_m+e_parallax_wd)

    return mask_parallax * mask_pm * mask_dist * mask_wd

def randomize(a, b, c, d, Ntot_wd):
    
    # Generate the permutation index array.
    # I'm randomizing the full catalog of white dwarfs.
    permutation = np.random.permutation(a.shape[0])
    
    # Shuffle the arrays by giving the permutation in the square brackets.
    shuffled_a = a[permutation]
    shuffled_b = b[permutation]
    shuffled_c = c[permutation]
    shuffled_d = d[permutation]
    
    #Return random proper motions for the number of white dwarfs I have
    return [shuffled_a[:Ntot_wd], shuffled_b[:Ntot_wd], 
            shuffled_c[:Ntot_wd], shuffled_d[:Ntot_wd]]

def calc_pca(params_m_0, params_wd_0, pwd, f_pwd, limits):
    '''
    Calculates the probability of chance alignment for the white dwarfs 
    m-dwarfs pairs by randomizing the proper motions of the white dwarfs. 
    Positions of both m and white dwarfs remain unchanged. 
    '''
        
    if(os.path.exists('Catalogs/wd_sources/gaiawd.fit')):
        w_dwarfs = Table.read('Catalogs/wd_sources/gaiawd.fit')
    else:
        url="ftp://cdsarc.u-strasbg.fr/pub/cats/J/MNRAS/482/4570/gaia2wd.dat.gz" 
        w_dwarfs = Table.read(url, readme="data/ReadMe", format="ascii.cds")
    
    Ntot_wd = len(pwd)
    N = 1000 #Number of repetition for the probability
    prob = np.zeros(Ntot_wd)
     
    #Set parameters of the white dwarfs original pairs
    parallax_wd_0, e_parallax_wd_0 = params_wd_0[0], params_wd_0[1]

    for i in range(N):
        #Randomize proper motions and errors using the full catalog of 
        #white dwarfs
        p = randomize(np.array(w_dwarfs['pmRA']), np.array(w_dwarfs['e_pmRA']), 
                      np.array(w_dwarfs['pmDE']), np.array(w_dwarfs['e_pmDE']),
                      Ntot_wd)
        
        pmra_wd_random,e_pmra_wd_random,pmdec_wd_random,e_pmdec_wd_random = p
        
        #Define new parameters list with the random proper motions
        params_wd_random = [parallax_wd_0, e_parallax_wd_0, 
                            pmra_wd_random, e_pmra_wd_random, 
                            pmdec_wd_random, e_pmdec_wd_random]
        
        #Decide if the pairs are co-movers or not with the random proper 
        #motions
        mask = get_mask_binaries(limits, params_m_0, params_wd_random, 
                                 pwd, f_pwd)

        #If there are any pairs, the probability of chance alignment increases 
        #for that m-dwarf
        prob[mask] += 1

    return np.array(prob)/N #Return normalized probability

def idenfy_good_wd_ages(m_w_co_movers):
    bp_rp = m_w_co_movers['BPmag']-m_w_co_movers['RPmag']
    g_abs = m_w_co_movers['Gmag'] - 5*(np.log10(1e3/m_w_co_movers['Plx'])-1)
    
    x=0.5
    model_da = np.loadtxt('Catalogs/Models/cooling_models/Table_Mass_'+str(x)+'_DA')
    g_model_da = model_da[:,34]
    bp_rp_model_da = model_da[:,35]-model_da[:,36]
    f_lim_up = interp1d(bp_rp_model_da,g_model_da)
    
    x=1.0
    model_da = np.loadtxt('Catalogs/Models/cooling_models/Table_Mass_'+str(x)+'_DA')
    g_model_da = model_da[:,34]
    bp_rp_model_da = model_da[:,35]-model_da[:,36]
    f_lim_down = interp1d(bp_rp_model_da,g_model_da)
    
    mask_discarded = np.logical_or(np.logical_or(
            f_lim_up(bp_rp) > g_abs,bp_rp > 0.9),
            f_lim_down(bp_rp) < g_abs)

    good_age = np.ones(len(bp_rp))
    good_age[mask_discarded] = 0
    
    return good_age
    