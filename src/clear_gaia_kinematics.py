#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from astropy.table import Table
from .astro import calc_number_single_stars

def mask_gaia_cuts(sample):
    '''
    Apply Gaia cuts to obtain the best kinematics. 
    The cuts defined in Kiman et al. 2019 were used.
    '''
    gaia_full = Table.read('Catalogs/literature_search_gaia_full.fits')
    phot_bp_rp_excess_factor = np.ones(len(sample))*np.nan
    visibility_periods_used = np.ones(len(sample))*np.nan
    astrometric_chi2_al = np.ones(len(sample))*np.nan
    astrometric_n_good_obs_al = np.ones(len(sample))*np.nan
    for i,x in enumerate(sample['source_id']):
        mask = x==gaia_full['source_id']
        if(len(gaia_full['phot_bp_rp_excess_factor'][mask])>0):
            phot_bp_rp_excess_factor[i] = gaia_full['phot_bp_rp_excess_factor'][mask][0]
            visibility_periods_used[i] = gaia_full['visibility_periods_used'][mask][0]
            astrometric_chi2_al[i] = gaia_full['astrometric_chi2_al'][mask][0]
            astrometric_n_good_obs_al[i] = gaia_full['astrometric_n_good_obs_al'][mask][0]

    sample['phot_bp_rp_excess_factor'] = phot_bp_rp_excess_factor
    sample['visibility_periods_used'] = visibility_periods_used
    sample['astrometric_chi2_al'] = astrometric_chi2_al
    sample['astrometric_n_good_obs_al'] = astrometric_n_good_obs_al
    
    uwe = np.sqrt(sample['astrometric_chi2_al']/(sample['astrometric_n_good_obs_al']-5))
    uwe_max = 1.2*np.array([np.max([1.4,np.exp(-0.2*(x-19.5))]) 
                            for x in sample['phot_g_mean_mag']])
    phot_bp_rp_excess_factor_max = 1.3 + 0.06*(sample['phot_bp_mean_mag']-sample['phot_rp_mean_mag'])**2
    mask_bp = (sample['phot_bp_mean_flux']/sample['phot_bp_mean_flux_error'])<=10
    phot_bp_rp_excess_factor_max[mask_bp]=99999
    
    #Gaia cuts
    mask_gaia = ((sample['parallax']/sample['parallax_error'] > 10)
             * (sample['visibility_periods_used'] > 8)
             * (uwe < uwe_max)
             * (sample['phot_bp_rp_excess_factor']<phot_bp_rp_excess_factor_max)
             * ((sample['phot_rp_mean_flux']/sample['phot_rp_mean_flux_error'])>10)
             * ((sample['phot_g_mean_flux']/sample['phot_g_mean_flux_error'])>10))
    
    #Register how many stars with good kinematics we have. We included a cut
    #to only have the number of single and compatible stars
    log_file = open('log.txt','a')
    mask_ha = ~np.isnan(sample['ewha'])
    n_kin = calc_number_single_stars(sample[mask_gaia*mask_ha])
    n_kin_tot = len(sample[mask_gaia*mask_ha])
    text = 'Literature search with good kinematics: {}\n'
    log_file.write(text.format(n_kin_tot))
    text = 'Number of single stars: {}\n'
    log_file.write(text.format(n_kin))
    log_file.flush()
    
    return sample[mask_gaia]