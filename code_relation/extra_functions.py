#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 16:34:55 2020

@author: rociokiman
"""

from astropy.table import Table
import numpy as np

def organize_table_format(ra,dec,source_id,ra_gaia,dec_gaia,pmra,pmra_error,
                          pmdec,pmdec_error,parallax,parallax_error,
                          phot_g_mean_mag,phot_rp_mean_mag,phot_bp_mean_mag,
                          ewha,ewha_error,ewha_all,ewha_error_all,
                          lhalbol,lhalbol_error,age,age_error_high,age_error_low,
                          group_num,group_name,source_num,source_ref):
    
    organized_table = Table()
    organized_table['ra'] = np.array(ra)
    organized_table['dec'] = np.array(dec)
    organized_table['gaia_source_id'] = np.array(source_id)
    organized_table['ra_gaia'] = np.array(ra_gaia)
    organized_table['dec_gaia'] = np.array(dec_gaia)
    organized_table['pmra']= np.array(pmra)
    organized_table['pmra_error']= np.array(pmra_error)
    organized_table['pmdec'] = np.array(pmdec)
    organized_table['pmdec_error'] = np.array(pmdec_error)
    organized_table['parallax'] = np.array(parallax)
    organized_table['parallax_error'] = np.array(parallax_error)
    organized_table['phot_g_mean_mag'] = np.array(phot_g_mean_mag)
    organized_table['phot_rp_mean_mag'] = np.array(phot_rp_mean_mag)
    organized_table['phot_bp_mean_mag'] = np.array(phot_bp_mean_mag)
    organized_table['ewha'] = np.array(ewha)
    organized_table['ewha_error'] = np.array(ewha_error)
    organized_table['ewha_all'] = np.array(ewha_all)
    organized_table['ewha_error_all'] = np.array(ewha_error_all)
    organized_table['lhalbol'] = np.array(lhalbol)
    organized_table['lhalbol_error'] = np.array(lhalbol_error)
    organized_table['age'] = np.array(age)
    organized_table['age_error_high'] = np.array(age_error_high)
    organized_table['age_error_low'] = np.array(age_error_low)
    organized_table['group_num'] = np.array(group_num)
    organized_table['group_name'] = np.array(group_name)
    organized_table['source_num'] = np.array(source_num)
    organized_table['source_ref'] = np.array(source_ref)
    
    return organized_table
    
    
    
    