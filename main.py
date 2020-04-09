#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.table import Table,vstack
import numpy as np
import code


###--------- M-dwarfs in moving groups ---------###

mg = Table.read('Catalogs/literature_search_mg.fits')

columns = [mg['ra'], mg['dec'], mg['source_id'], mg['ra_x'], mg['dec_x'],
           mg['pmra'], mg['pmra_error'], mg['pmdec'], mg['pmdec_error'],
           mg['parallax'], mg['parallax_error'], mg['phot_g_mean_flux'],
           mg['phot_g_mean_flux_error'], mg['phot_g_mean_mag'],
           mg['phot_rp_mean_flux'], mg['phot_rp_mean_flux_error'],
           mg['phot_rp_mean_mag'], mg['phot_bp_mean_flux'],
           mg['phot_bp_mean_flux_error'], mg['phot_bp_mean_mag'], mg['ewha'],
           mg['ewha_error'], mg['ewha_all'], mg['ewha_error_all'], 
           mg['lhalbol'], mg['lhalbol_error'], mg['age']*1e6, 
           (mg['age_error']*1e6)/2, (mg['age_error']*1e6)/2, mg['group_num'],
           mg['group_name'], mg['source_num'], mg['source_ref']]

m_dwarfs_mg = code.organize_table_format(columns)

del mg

###--------- M-dwarfs co-moving with a White dwarf ---------###

#Load sample of m-dwarfs
m_dwarfs_not_mg = Table.read('Catalogs/literature_search_not_mg.fits')

#Cross-mach m-dwarf sample with white dwarf sample to find all the pairs in
#a 10 arcmin radius
m_dwarfs_pairs,w_dwarfs_pairs,_ = code.cross_match_to_white_dwarfs(m_dwarfs_not_mg)

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
mask_comovers_all = code.get_mask_binaries(limits,params_m,params_wd,
                                            w_dwarfs_pairs['Pwd'],
                                            w_dwarfs_pairs['f_Pwd'])

#Calculate probability of chance alignment for each pair
prob_chance_align = code.calc_pca(params_m, params_wd,
                                   w_dwarfs_pairs['Pwd'], 
                                   w_dwarfs_pairs['f_Pwd'],
                                   limits)

#Final mask for pairs includes selected co-movers which have a low probability 
#of chance alignment
mask_comovers = mask_comovers_all * (prob_chance_align <= 0.01)

info = 'Number of m-dwarfs white dwarfs pairs: {}'
N = len(np.array(m_dwarfs_pairs['parallax'])[mask_comovers_all])
print(info.format(N))
info = 'Number of m-dwarfs with prob of chance alignment > 0.01: {}'
mask_info = mask_comovers_all * (prob_chance_align > 0.01)
N = len(np.array(m_dwarfs_pairs['parallax'])[mask_info])
print(info.format(N))
info = 'Number of m-dwarfs white dwarfs pairs with prob ' \
       'of chance alignment <= 0.01: {}'
N = len(np.array(m_dwarfs_pairs['parallax'])[mask_comovers])
print(info.format(N))

#Remove nans to calculate ages in the next step.
mask_nan_teff_logg = ~np.isnan(w_dwarfs_pairs['TeffH']+w_dwarfs_pairs['loggH'])

info = 'Number of m-dwarfs white dwarfs pairs I will calculate age for: {}'
N = len(np.array(m_dwarfs_pairs['parallax'])[mask_comovers*mask_nan_teff_logg])
print(info.format(N))

#Define the new sample of m-dwarfs that have a white dwarf co-moving
m_co_movers = m_dwarfs_pairs[mask_comovers*mask_nan_teff_logg]
w_co_movers = w_dwarfs_pairs[mask_comovers*mask_nan_teff_logg]

w_co_movers.write('Catalogs/wdm_binaries.fits', format = 'fits', 
                  overwrite = True)

#Calculate total age, cooling age, main sequence age, initial mass 
#and final mass for the white dwarfs

result_w_ages = code.calc_ages_wdm_binaries(w_co_movers)

#Save table of m dwarfs co-moving with white dwarfs in a nice format
N_final = len(m_co_movers)

columns = [m_co_movers['ra'], m_co_movers['dec'], m_co_movers['source_id'],
           m_co_movers['ra_x'], m_co_movers['dec_x'], m_co_movers['pmra'], 
           m_co_movers['pmra_error'], m_co_movers['pmdec'],
           m_co_movers['pmdec_error'], m_co_movers['parallax'],
           m_co_movers['parallax_error'], m_co_movers['phot_g_mean_flux'],
           m_co_movers['phot_g_mean_flux_error'], 
           m_co_movers['phot_g_mean_mag'], m_co_movers['phot_rp_mean_flux'],
           m_co_movers['phot_rp_mean_flux_error'], 
           m_co_movers['phot_rp_mean_mag'], m_co_movers['phot_bp_mean_flux'],
           m_co_movers['phot_bp_mean_flux_error'], 
           m_co_movers['phot_bp_mean_mag'], m_co_movers['ewha'], 
           m_co_movers['ewha_error'], m_co_movers['ewha_all'], 
           m_co_movers['ewha_error_all'], m_co_movers['lhalbol'],
           m_co_movers['lhalbol_error'], result_w_ages[:,6],
           result_w_ages[:,7], result_w_ages[:,8],
           np.zeros(N_final), np.array(['WD' for i in range(N_final)]),
           m_co_movers['source_num'], m_co_movers['source_ref']]

m_dwarfs_wd = code.organize_table_format(columns)

###--------- Final Age-calibrators sample ---------###

#Combine results 
age_calibrators = vstack([m_dwarfs_mg,m_dwarfs_wd])

print('Total number of age calibrators: {}'.format(len(age_calibrators)))

###--------- Calculate LHalphaLbol ---------###

lhalbol,lhalbol_error = code.calc_lhalbol(age_calibrators['ewha'],
                                           age_calibrators['ewha_error'],
                                           age_calibrators['phot_g_mean_mag']-
                                           age_calibrators['phot_rp_mean_mag'])

age_calibrators['lhalbol'] = lhalbol
age_calibrators['lhalbol_error'] = lhalbol_error

age_calibrators.write('Catalogs/age_calibrators_bayes.fits', overwrite=True)




