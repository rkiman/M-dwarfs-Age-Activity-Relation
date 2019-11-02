#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 15:32:02 2019

@author: rociokiman
"""
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np

def hstodeg(ra_hs,dec_hs):
    ra,dec = [],[]
    for i in range(len(ra_hs)):
        x,y = ra_hs[i],dec_hs[i]
        if(str(x)=='--' or str(y)=='--'):
            ra.append(np.nan)
            dec.append(np.nan)        
        else:
            radec_hs_i = x + ' ' +  y
            c = SkyCoord(radec_hs_i, unit=(u.hourangle, u.deg))
            ra.append(c.ra.deg)
            dec.append(c.dec.deg)
    return np.array(ra),np.array(dec)