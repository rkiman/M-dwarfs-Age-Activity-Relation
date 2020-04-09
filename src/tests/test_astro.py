#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from src.astro import WISE_id_to_ra_dec
import numpy as np

def test_WISE_id_to_ra_dec():
    ra,dec = WISE_id_to_ra_dec(np.array(['J081417.67+025319.4',
                                         'J082543.32-002911.4']))
    assert np.isclose(ra[0],123.5736250) 
    assert np.isclose(dec[0],2.8887222)
    assert np.isclose(ra[1],126.4305000)       
    assert np.isclose(dec[1],-0.4865000)

