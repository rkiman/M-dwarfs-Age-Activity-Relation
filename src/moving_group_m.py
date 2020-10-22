#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import sys
bpath='/Users/rociokiman/Documents/Gaia-Cupid/ActivityAgeRelation/banyan_sigma'
sys.path.append(bpath)
from banyan_sigma import banyan_sigma
from scipy.io.idl import readsav
from astropy.table import Table
from .astro import organize_table_format,calc_number_single_stars

def compile_m_moving_groups_sample(ls_compatible):
    '''
    From the literature search already checked that it is compatible with the
    Kiman et al. 2019 catalog, compiles a table with moving group members from
    than list using Banyan Sigma (Gagné,J. et al., Astrophys.J. 856, 23(2018)).
    '''
    
    #Load reference of known moving groups in Banyan
    mg_ref = Table.read('data/moving_groups_ref.csv',format='csv')
    
    #columns needed in this part of the code
    ra = np.array([float(x) for x in ls_compatible['ra_x']])
    dec = np.array([float(x) for x in ls_compatible['dec_x']])
    pmra = np.array([float(x) for x in ls_compatible['pmra']])
    pmra_error = np.array([float(x) for x in ls_compatible['pmra_error']])
    pmdec = np.array([float(x) for x in ls_compatible['pmdec']])
    pmdec_error = np.array([float(x) for x in ls_compatible['pmdec_error']])
    parallax = np.array([float(x) for x in ls_compatible['parallax']])
    parallax_error = [float(x) for x in ls_compatible['parallax_error']]
    parallax_error = np.array(parallax_error)
    rv = np.array([float(x) for x in ls_compatible['radial_velocity']])
    rv_error = [float(x) for x in ls_compatible['radial_velocity_error']]
    rv_error = np.array(rv_error)
    group_name = np.array([str(x) for x in ls_compatible['group_name']])
    source_id = np.array([str(x) for x in ls_compatible['source_id']])
    
    #Define mask to run banyan correctly
    mask_run_banyan = (~np.isnan(ra+dec+pmra+pmra_error+pmdec+pmdec_error+
                                 parallax+parallax_error))

    #Run banyan
    print('Running Banyan')
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
    best_hyp = np.array(result['BEST_HYP']).reshape(len(result['BEST_HYP']),)

    #Get Praesepe members which is not included in banyan yet.
    #Jonathan Gagne run Banyan on this sample for me and here are the results
    path = 'Catalogs/Sources/praesepe_sample.sav'
    #extracting best young asociationg and probability of being member
    results_pra = readsav(path)
    ya_prob_pra = results_pra['out']['YA_PROB']
    best_ya_pra = results_pra['out']['BEST_YA']
    best_hyp_pra = results_pra['out']['BEST_HYP']
    source_id_pra = results_pra['input']['source_id']
    
    #Columns needed for only the good kinematic stars
    group_run_banyan = group_name[mask_run_banyan]
    source_id_run_banyan = source_id[mask_run_banyan]
    
    #Adjust PRA members with the results from Jonathan's run from Banyan in
    #PRA
    print('Adjusting PRA members')
    prob_ya,best_ya,best_hyp = adjust_pra_members(group_run_banyan,
                                                  source_id_run_banyan,prob_ya,
                                                  best_ya,best_hyp,
                                                  source_id_pra,ya_prob_pra,
                                                  best_ya_pra,best_hyp_pra)
    
    best_ya = np.array([str(x) for x in best_ya])
    best_hyp = np.array([str(x) for x in best_hyp])
    
    #Mask for high likelihood members according to banyan
    print('Making mask for members')
    mask_and_nums = get_membership_mask(group_run_banyan,source_id_run_banyan,
                                        result,prob_ya,best_ya,best_hyp)
    
    mask_membership,mem_low,field,uvw_sep_high,good_mem,arg = mask_and_nums
    
    #Save fits file with all the stars and their group classification to 
    #identify new members.
    all_groups = ls_compatible[mask_run_banyan]
    all_groups['prob_ya'] = prob_ya
    all_groups['best_ya'] = best_ya
    all_groups['best_hyp'] = best_hyp
    all_groups['mem_low'] = mem_low
    all_groups['field'] = field
    all_groups['uvw_sep_high'] = uvw_sep_high
    all_groups['good_mem'] = good_mem
    all_groups['arg'] = arg
    
    all_groups.write('Catalogs/literature_search_all_groups.fits',
                     format='fits',overwrite=True)

    #Register number of stars identified as members, not members and why
    mask_ha = ~np.isnan(all_groups['ewha'])
    subsample = all_groups[(all_groups['mem_low']==1)*mask_ha]
    n_mem_low = calc_number_single_stars(subsample)
    subsample = all_groups[(all_groups['field']==1)*mask_ha]
    n_field = calc_number_single_stars(subsample)
    subsample = all_groups[(all_groups['uvw_sep_high']==1)*mask_ha]
    n_uvw_sep_high = calc_number_single_stars(subsample)
    subsample = all_groups[(all_groups['good_mem']==1)*mask_ha]
    n_good_mem_tot = len(subsample)
    n_good_mem = calc_number_single_stars(subsample)
    subsample = all_groups[(all_groups['arg']==1)*mask_ha]
    n_arg = calc_number_single_stars(subsample)

    log_file = open('log.txt','a')
    log_file.write('Number of single compatible stars not in moving group because the prob of the group was too low:{}\n'.format(n_mem_low))
    log_file.write('Number of single compatible stars not in moving group because they were classified as field stars:{}\n'.format(n_field))
    log_file.write('Number of single compatible stars not in moving group because they had a distance in UVW > 5km/s:{}\n'.format(n_uvw_sep_high))
    log_file.write('Number of single compatible stars not in moving group because they were part of argus:{}\n'.format(n_arg))
    log_file.write('Number of high likelihood members:{}\n'.format(n_good_mem_tot))
    log_file.write('Number of single high likelihood members:{}\n'.format(n_good_mem))
    log_file.flush()
    
    #Create table with the true members accorsing to Banyan
    mg_sample = ls_compatible[mask_run_banyan][mask_membership]
    
    mg_sample['best_ya'] = np.array([str(x) for x in best_ya[mask_membership]])
    
    #Create sample of stars that don't belong to a moving group:
    m_dwarfs_not_mg = ls_compatible[mask_run_banyan][~mask_membership]
    
    #Correcting groups that don't agree with banyan
    for i in range(len(mg_sample)):
        if(mg_sample['group_name'][i]!=mg_sample['best_ya'][i]):
            mask = mg_sample['best_ya'][i] == mg_ref['name']
            mg_sample['group_name'][i] = mg_ref['name'][mask][0]
            mg_sample['group_num'][i] = mg_ref['group_num'][mask][0]
            mg_sample['age'][i] = mg_ref['age'][mask][0]
            mg_sample['age_error'][i] = mg_ref['age_error'][mask][0]

    #Organize table to use in future steps
    columns = [mg_sample['ra'],mg_sample['dec'],mg_sample['spt'],
               mg_sample['source_id'], 
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
               mg_sample['star_index'],
               mg_sample['source_num'], mg_sample['source_ref']]

    m_dwarfs_mg = organize_table_format(columns)
    
    print('Saving literature_search_mg.fits')
    m_dwarfs_mg.write('Catalogs/literature_search_mg.fits',format='fits',
                      overwrite=True)
    
    return m_dwarfs_mg,m_dwarfs_not_mg

def remove(string): 
    return string.replace(" ", "")

def adjust_pra_members(group,source_id,prob_ya,best_ya,best_hyp,
                       source_id_pra,ya_prob_pra,best_ya_pra,best_hyp_pra):
    for i in range(len(prob_ya)):
        if(group[i]=='PRA'):
            mask_pra_i = [int(source_id[i])==int(x) for x in source_id_pra]
            mask_pra_i = np.array(mask_pra_i)
            x = ya_prob_pra[mask_pra_i][0]
            y = best_ya_pra[mask_pra_i][0]
            z = best_hyp_pra[mask_pra_i][0]
            #find all the stars with the same gaia id 
            mask_source_id = [int(source_id[i])==int(x) for x in source_id]
            mask_source_id = np.array(mask_source_id)
            if(np.isnan(x)):
                prob_ya[i] = np.nan
                best_ya[i] = np.nan
                best_hyp[i] = np.nan
            else:
                prob_ya[mask_source_id] = x
                best_ya[mask_source_id] = y.decode("utf-8") 
                best_hyp[mask_source_id] = z.decode("utf-8") 
                if(best_ya[i] == 'PRAE'):
                    best_ya[mask_source_id] = 'PRA'
                if(best_hyp[i] == 'PRAE'):
                    best_hyp[mask_source_id] = 'PRA'
    return prob_ya,best_ya,best_hyp

def get_membership_mask(group,source_id,result,prob_ya,best_ya,best_hyp):
    n_all = len(prob_ya)
    n_mem_low = np.zeros(n_all)
    n_field = np.zeros(n_all)
    n_uvw_sep_high = np.zeros(n_all)
    n_good_mem = np.zeros(n_all)
    n_arg = np.zeros(n_all)
    highprob = []
    for i in range(len(prob_ya)):   
        x = prob_ya[i]
        y = best_ya[i]
        z = best_hyp[i]
        if(x > 0.9):
            error_msg = 'Check why best_ya: {}!= best_hyp {}'.format(y,z)
            assert y==z, error_msg
            if(y!='PRA' and result[str.encode(z)]['UVW_SEP'][i]>5):
                n_uvw_sep_high[i] = 1
                highprob.append(False)
            elif(y=='ARG'):
                n_arg[i] = 1
                highprob.append(False)
            else:
                highprob.append(True)
                n_good_mem[i] = 1
        elif(x <= 0.9):
            highprob.append(False)
            if(best_hyp[i]=='FIELD'):
                n_field[i] = 1
            elif(best_hyp[i]!='FIELD'):
                n_mem_low[i] = 1
    return np.array(highprob),n_mem_low,n_field,n_uvw_sep_high,n_good_mem,n_arg