#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from src.astro import WISE_id_to_ra_dec,calc_number_single_stars
import numpy as np
from astropy.table import Table

def test_WISE_id_to_ra_dec():
    ra,dec = WISE_id_to_ra_dec(np.array(['J081417.67+025319.4',
                                         'J082543.32-002911.4']))
    assert np.isclose(ra[0],123.5736250) 
    assert np.isclose(dec[0],2.8887222)
    assert np.isclose(ra[1],126.4305000)       
    assert np.isclose(dec[1],-0.4865000)
    
def test_calc_number_single_stars():
    test_table = Table()
    same_star_single = np.array([1, 2, 3, 4, 5, 6])
    n_nan = 8
    n_single_true = len(same_star_single) + n_nan
    same_star = []
    for x in same_star_single:
        n = np.random.randint(1,5)
        for i in range(n):
            same_star.append(x)
    for i in range(n_nan):
        same_star.append(np.nan)
    test_table['same_star'] = np.array(same_star)
    n_single = calc_number_single_stars(test_table)
    
    assert n_single == n_single_true

