#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from src.cross_validation import calc_cross_validation
import numpy as np


def f(x,a,b):
    return a*x+b

def g(x,*p):
    return np.polyval(p,x)

def test_calc_cross_validation():
    a,b = 1,2
    n = 10
    sigma = 1
    x = np.linspace(0,10,n)
    y = np.random.normal(f(x,a,b),sigma)
    y_err = np.ones(n)*sigma
    ini_params1 = np.array([1,1])
    score1 = calc_cross_validation(x,y,y_err,f,ini_params1)
    ini_params2 = np.array([1,1,1,1,1])
    score2 = calc_cross_validation(x,y,y_err,g,ini_params2)
    
    assert score1 < score2
    

