#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.table import Table
import numpy as np
import pkg_resources

def organize_banyan_sample():
    
    #Load young moving groups catalog from banyan sigma
    path = '../Catalogs/YMG_Bona_Fide.csv'
    filepath = pkg_resources.resource_filename(__name__, path)
    
    mg_catalog = Table.read(filepath)
    
    print(mg_catalog)
#def confirm_active_m_dwarfs(cluster=''):
#    if(cluster=='HYA'):
#        table_hya = 
        
        
        
organize_banyan_sample()

