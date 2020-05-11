#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import sys
sys.path.append('/Users/rociokiman/Documents/Gaia-Cupid/ActivityAgeRelation/banyan_sigma')
from banyan_sigma import banyan_sigma
from scipy.io.idl import readsav
from astropy.table import Table
import os
from .astro import organize_table_format

def compile_m_moving_groups_sample(ls_compatible):
    '''
    From the literature search already checked that it is compatible with the
    Kiman et al. 2019 catalog, compiles a table with moving group members from
    than list using Banyan Sigma (Gagné,J. et al., Astrophys.J. 856, 23(2018)).
    '''
    
    #Load reference of known moving groups in Banyan
    mg_ref = Table.read('data/moving_groups_ref.csv',format='csv')
    
    N = len(ls_compatible)
    
    ra = np.array([float(x) for x in ls_compatible['ra_x']])
    dec = np.array([float(x) for x in ls_compatible['dec_x']])
    pmra = np.array([float(x) for x in ls_compatible['pmra']])
    pmra_error = np.array([float(x) for x in ls_compatible['pmra_error']])
    pmdec = np.array([float(x) for x in ls_compatible['pmdec']])
    pmdec_error = np.array([float(x) for x in ls_compatible['pmdec_error']])
    parallax = np.array([float(x) for x in ls_compatible['parallax']])
    parallax_error = np.array([float(x) for x in ls_compatible['parallax_error']])
    rv = np.array([float(x) for x in ls_compatible['radial_velocity']])
    rv_error = np.array([float(x) for x in ls_compatible['radial_velocity_error']])

    
    #Define mask to run banyan correctly
    mask_run_banyan = (~np.isnan(ra+dec+pmra+pmra_error+pmdec+pmdec_error+
                                 parallax+parallax_error)
                       * (parallax/parallax_error > 8))
    
    n_kin = len(ra[mask_run_banyan])
    print('Number of not accretors with good kinematics: {}'.format(n_kin))
    
    #Run banyan
    result = banyan_sigma(ra=ra[mask_run_banyan],dec=dec[mask_run_banyan], 
                          pmra=pmra[mask_run_banyan],
                          pmdec=pmdec[mask_run_banyan], 
                          epmra=pmra_error[mask_run_banyan], 
                          epmdec=pmdec_error[mask_run_banyan],
                          rv=rv[mask_run_banyan],erv=rv_error[mask_run_banyan], 
                          plx=parallax[mask_run_banyan], 
                          eplx=parallax_error[mask_run_banyan])
    #Results from banyan
    prob_ya = np.array(result['YA_PROB']).reshape(len(result['YA_PROB']),)
    best_ya = np.array(result['BEST_YA']).reshape(len(result['BEST_YA']),)
    
    #for i in range(len(prob_ya)):
    #    if(str(best_ya[i])=='ARG'):
    #        prob_ya[i] = 0
            
    #Mask for high likelihood members accordin to banyan
    highprob = prob_ya > 0.9
    
    #Get Praesepe members which is not included in banyan yet.
    #Jonathan Gagne run Banyan on this sample for me and here are the results
    path = 'Catalogs/rocio_praesepe_sample.sav'
    #extracting best young asociationg and probability of being member
    results_pra = readsav(path,verbose=1)
    #best_ya_bytes = results['out']['best_ya']
    #best_ya_pra = np.array([x.decode("utf-8") for x in best_ya_bytes])
    ya_prob_pra = results_pra['out']['YA_PROB']
    source_id_pra = results_pra['input']['source_id'][ya_prob_pra>0.9]
    
    #Mask for PRA members
    bf_pra = np.array([True if x in source_id_pra else False 
                       for x in ls_compatible['source_id']])
    
    #Ajust prob for Praesepe members
    for x,y in zip(source_id_pra,ya_prob_pra[ya_prob_pra>0.9]):
        mask = ls_compatible['source_id'][mask_run_banyan] == x
        prob_ya[mask] = y
    
    #Making sure the true members in Praesepe have the right age
    for i in range(N):
        if(remove(ls_compatible['group_name'][i])=='PRA' and bf_pra[i]):
            mask_ref = mg_ref['name']=='PRA'
            ls_compatible['group_name'][i] = 'PRA'
            ls_compatible['group_num'][i] = mg_ref['group_num'][mask_ref][0]
            ls_compatible['age'][i] = mg_ref['age'][mask_ref][0]
            ls_compatible['age_error'][i] = mg_ref['age_error'][mask_ref][0]
        
            
    #Define mask for high likelihood members including praesepe        
    mask_membership = np.logical_or(highprob,bf_pra[mask_run_banyan])
    
    mg_sample = ls_compatible[mask_run_banyan][mask_membership]
    mg_sample['ya_prob'] = prob_ya[mask_membership]
    mg_sample['best_ya'] = np.array([str(x) for x in best_ya[mask_membership]])
    
    n_mg = len(mg_sample[~np.isnan(mg_sample['ewha'])])
    print("Number of high like mem: {}".format(n_mg))
    
    #Create sample of stars that don't belong to a moving group:
    m_dwarfs_not_mg = ls_compatible[mask_run_banyan][~mask_membership]
    if(~os.path.exists('Catalogs/literature_search_not_mg.fits')):
        print('Saving literature_search_not_mg.fits')
        m_dwarfs_not_mg.write('Catalogs/literature_search_not_mg.fits',
                              format='fits')
    
    #Correcting groups that don't agree with banyan
    for i in range(len(mg_sample)):
        if(mg_sample['group_name'][i]!=mg_sample['best_ya'][i]):
            if(mg_sample['group_name'][i]!='PRA'):
                mask = mg_sample['best_ya'][i] == mg_ref['name']
                mg_sample['group_name'][i] = mg_ref['name'][mask][0]
                mg_sample['group_num'][i] = mg_ref['group_num'][mask][0]
                mg_sample['age'][i] = mg_ref['age'][mask][0]
                mg_sample['age_error'][i] = mg_ref['age_error'][mask][0]
            elif(mg_sample['group_name'][i]=='PRA'):
                mg_sample['best_ya'][i] = 'PRA'
                mg_sample['ya_prob'][i] = np.nan

    #Organize table to use in future steps
    columns = [mg_sample['ra'],mg_sample['dec'],mg_sample['source_id'], 
               mg_sample['ra_x'],mg_sample['dec_x'],mg_sample['pmra'],
               mg_sample['pmra_error'],mg_sample['pmdec'],
               mg_sample['pmdec_error'],mg_sample['parallax'], 
               mg_sample['parallax_error'],mg_sample['phot_g_mean_flux'],
               mg_sample['phot_g_mean_flux_error'], 
               mg_sample['phot_g_mean_mag'],mg_sample['phot_rp_mean_flux'], 
               mg_sample['phot_rp_mean_flux_error'],
               mg_sample['phot_rp_mean_mag'],mg_sample['phot_bp_mean_flux'],
               mg_sample['phot_bp_mean_flux_error'],
               mg_sample['phot_bp_mean_mag'],
               mg_sample['g_corr'],mg_sample['rp_corr'],
               mg_sample['ewha'],
               mg_sample['ewha_error'],mg_sample['ewha_all'], 
               mg_sample['ewha_error_all'],mg_sample['lhalbol'], 
               mg_sample['lhalbol_error'], mg_sample['age']*1e6,
               (mg_sample['age_error']*1e6)/2,(mg_sample['age_error']*1e6)/2, 
               mg_sample['group_num'],mg_sample['group_name'], 
               mg_sample['source_num'], mg_sample['source_ref']]

    m_dwarfs_mg = organize_table_format(columns)
    if(~os.path.exists('Catalogs/literature_search_mg.fits')):
        print('Saving literature_search_mg.fits')
        m_dwarfs_mg.write('Catalogs/literature_search_mg.fits',format='fits')
    
    return m_dwarfs_mg,m_dwarfs_not_mg

def remove(string): 
    return string.replace(" ", "")

