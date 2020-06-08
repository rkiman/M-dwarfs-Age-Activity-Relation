#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.table import Table
import src 

#binaries = Table.read('Catalogs/wdm_binaries.fits')
#src.make_table_for_wd(binaries)
#src.make_table_wd_ages(binaries)
    
age_calibrators = Table.read('Catalogs/age_calibrators_bayes.fits')
src.make_table_summary_age_calibrators(age_calibrators)
    
#accretors = Table.read('Catalogs/literature_search_accretors.fits')
#src.make_table_accretors(accretors)

