#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.table import Table
import numpy as np
import pkg_resources
from astropy.coordinates import SkyCoord
from astropy import units as u
#from .astro import WISE_id_to_ra_dec
import sys
sys.path.append('/Users/rociokiman/Documents/Gaia-Cupid/ActivityAgeRelation/banyan_sigma')
from banyan_sigma import banyan_sigma
from scipy.io.idl import readsav


def str_to_float(array):
    n = len(array)
    array_float = np.ones(n)*np.nan
    for i in range(n):
        if(array[i]!='...'):
           array_float[i] = float(array[i])
    return array_float

def organize_banyan_sample():
    
    #Load young moving groups catalog from banyan sigma
    path = '../Catalogs/YMG_Bona_Fide.csv'
    filepath = pkg_resources.resource_filename(__name__, path)
    
    mg_catalog = Table.read(filepath)

    #wise_name = np.array(mg_catalog['col17'][2:])
    #ra,dec = WISE_id_to_ra_dec(wise_name)
    ra = np.array([float(x) for x in mg_catalog['col62'][2:]])
    dec = np.array([float(x) for x in mg_catalog['col63'][2:]])
    
    organized_banyan = Table()

    organized_banyan['ra'] = ra
    organized_banyan['dec'] = dec    
    organized_banyan['pmra'] = str_to_float(mg_catalog['col78'][2:])
    organized_banyan['e_pmra'] = str_to_float(mg_catalog['col79'][2:])
    organized_banyan['pmdec'] = str_to_float(mg_catalog['col80'][2:])
    organized_banyan['e_pmdec'] = str_to_float(mg_catalog['col81'][2:])
    organized_banyan['rv'] = str_to_float(mg_catalog['col90'][2:])
    organized_banyan['e_rv'] = str_to_float(mg_catalog['col91'][2:])
    organized_banyan['parallax'] = str_to_float(mg_catalog['col95'][2:])
    organized_banyan['e_parallax'] = str_to_float(mg_catalog['col96'][2:])
    organized_banyan['group'] = np.array(mg_catalog[';papers_to_add'][2:])
    organized_banyan['mem_type'] = np.array(mg_catalog['TODO'][2:])

    return organized_banyan

def crossmatch_sample_w_banyan(ra,dec):
    organized_banyan = organize_banyan_sample()
    
    sample = SkyCoord(ra=np.array(ra)*u.degree, 
                      dec=np.array(dec)*u.degree)
    banyan = SkyCoord(ra=np.array(organized_banyan['ra'])*u.degree,
                      dec=np.array(organized_banyan['dec'])*u.degree)
    
    #Cross-match the two samples with a 5arcsec radius
    idx_sample, idx_banyan, d2d,_ = banyan.search_around_sky(sample,5*u.arcsec)
    
    #Transform separation to arcmin to check in the test file
    separation = d2d.to(u.arcsec)/u.arcsec
    
    return idx_sample, organized_banyan[idx_banyan], separation

def confirm_active_m_dwarfs(cluster=''):
    
    if(cluster=='HYA'):
        path = '../Catalogs/Sources/Douglas2014_hyades_gaia.csv'
        filepath = pkg_resources.resource_filename(__name__, path)
        douglas2014 = Table.read(filepath)        
        ra = np.array(douglas2014['ra'])
        dec = np.array(douglas2014['dec'])
        parallax = np.array(douglas2014['parallax'])
        parallax_error = np.array(douglas2014['parallax_error'])
        pmra = np.array(douglas2014['pmra'])
        pmra_error = np.array(douglas2014['pmra_error'])
        pmdec = np.array(douglas2014['pmdec'])
        pmdec_error = np.array(douglas2014['pmdec_error'])

        N = len(ra)

        mask_run_banyan = (~np.isnan(ra+dec+pmra+pmra_error+pmdec+pmdec_error+
                                     parallax+parallax_error)
                          * (parallax/parallax_error > 8))
            
        result = banyan_sigma(ra=ra[mask_run_banyan], dec=dec[mask_run_banyan],
                              pmra=pmra[mask_run_banyan],
                              pmdec=pmdec[mask_run_banyan], 
                              epmra=pmra_error[mask_run_banyan],
                              epmdec=pmdec_error[mask_run_banyan],
                              plx=parallax[mask_run_banyan],
                              eplx=parallax_error[mask_run_banyan])
    
        prob_ya = np.array(result['YA_PROB']).reshape(len(result['YA_PROB']),)
        best_ya = np.array(result['BEST_YA']).reshape(len(result['BEST_YA']),)
        
        prob_mg_mem = np.ones(N)*np.nan
        best_ya_mem = np.array(['none' for x in range(N)])
        
        prob_mg_mem[mask_run_banyan] = prob_ya
        best_ya_mem[mask_run_banyan] = best_ya
        
        douglas2014['prob_ya'] = prob_mg_mem
        douglas2014['best_ya'] = best_ya_mem
        
        return douglas2014
    
    if(cluster=='PRA'):
        path = '../Catalogs/Sources/Douglas2014_praesepe_gaia.csv'
        filepath = pkg_resources.resource_filename(__name__, path)
        douglas2014 = Table.read(filepath)
        
        #Jonathan Gagne run Banyan on this sample for me and here are the
        #results
        path = '../Catalogs/praesepe_sample_banyan_sigma.sav'
        filepath = pkg_resources.resource_filename(__name__, path)
        
        #extracting best young asociationg and probability of being member
        results = readsav(filepath,verbose=1)
        best_ya_bytes = results['out']['best_ya']
        best_ya = np.array([x.decode("utf-8") for x in best_ya_bytes])
        ya_prob = results['out']['YA_PROB']
        
        douglas2014['best_ya'] = best_ya
        douglas2014['ya_prob'] = ya_prob
        
        return douglas2014
        
    
def confirm_active_m_dwarfs1(cluster=''):
    
    if(cluster=='HYA'):
        path = '../Catalogs/Sources/Douglas2014yades.fit'
        filepath = pkg_resources.resource_filename(__name__, path)
        
        table_hya = Table.read(filepath)

        results = crossmatch_sample_w_banyan(np.array(table_hya['_RAJ2000']),
                                             np.array(table_hya['_DEJ2000']))
        
        idx_sample, banyan, separation = results
        
        table_hya = table_hya[idx_sample]
        
        for x in banyan.colnames:
            table_hya[x] = banyan[x]
        
        mem_type_num = np.ones(len(table_hya))*np.nan
        
        for i,mtype in zip(range(5),['R','LM','CM','HM','BF']):
            mask = table_hya['mem_type'] == mtype
            mem_type_num[mask] = i
        
        table_hya['mem_type_num'] = mem_type_num
        table_hya['separation'] = separation
            
        return table_hya

