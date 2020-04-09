#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from src.WhiteDwarfsComovers import cross_match_to_white_dwarfs,calc_pca
from astropy.table import Table
import numpy as np


def test_cross_match_to_white_dwarfs():
    
    m_dwarf = Table()
    m_dwarf['ra'] = np.array([2.8825499999999997,2.9710166666666664,
           2.9851666666666663,4.0608958333333325])
    m_dwarf['dec'] = np.array([59.14442777777777,22.984545833333332,33.05468,
           19.86043611111111])
    
    _,_,separation = cross_match_to_white_dwarfs(m_dwarf)
    
    assert all(separation <= 10) #only matches with a 10arcmin separation
    
def test_calc_pca():
    
    params_m_0 = [np.array([14.3647,5.1969]),np.array([0.0404,0.0832]),
                  np.array([-20.777,-36.256]),np.array([0.06,0.118]),
                  np.array([-106.458,-13.225]),np.array([0.072,0.072])]
                  
    params_wd_0 = [np.array([14.515, 5.115]),np.array([0.1599,0.2005]),
                   np.array([-19.859,-35.33 ]),np.array([0.267,0.304]),
                   np.array([-107.222,-12.697]),np.array([0.33,0.188])]
    
    pwd, f_pwd = np.array([0.99577,0.99896]), np.array([0,0])
    
    limits = [8, #parallax snr cut for m dwarfs
              4, #parallax snr cut for white dwarfs
              3, #sigma difference for pmra and pmdec
              3] #sigma difference for parallax
    
    prob_chance_align = calc_pca(params_m_0, params_wd_0, pwd, f_pwd, limits)
    
    assert all(prob_chance_align <= 0.01)
    
    