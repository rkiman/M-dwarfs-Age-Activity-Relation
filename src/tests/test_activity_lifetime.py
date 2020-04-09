#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.table import Table
import numpy as np
import pkg_resources
from scipy.io.idl import readsav


def text_confirm_active_m_dwarfs():
    path = '../Catalogs/Sources/Douglas2014_praesepe_gaia.csv'
    filepath = pkg_resources.resource_filename(__name__, path)
    douglas2014 = Table.read(filepath)
    
    #Jonathan Gagne run Banyan on this sample for me and here are the
    #results
    path = '../Catalogs/praesepe_sample_banyan_sigma.sav'
    filepath = pkg_resources.resource_filename(__name__, path)
    results = readsav(filepath,verbose=1)
    
    i = int(np.random.uniform(0,len(douglas2014[1].data['ra']),10))
    
    assert(all(results['input']['ra'][i] == douglas2014['ra'][i]))
    assert(all(results['input']['dec'][i] == douglas2014['dec'][i]))
