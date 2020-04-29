#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import curve_fit

def calc_cross_validation(x,y,y_err,f,ini_params):
    score = 0
    N = len(y)
    
    for k in range(N):
        x_k = x[k]
        y_k = y[k]
        y_err_k = y_err[k]
        x_train = np.delete(x,k)
        y_train = np.delete(y,k)
        y_err_train = np.delete(y_err,k)
        
        popt,pcov=curve_fit(f,x_train,y_train,p0=ini_params,sigma=y_err_train)
        y_pred = f(x_k,*popt)

        score += calc_cross_validation_score(y_pred,y_k,y_err_k)
    return score
        
def calc_cross_validation_score(y_pred,y_k,y_err_k):
    return np.sum((y_k-y_pred)**2 / y_err_k**2)
    
        
'''
x = np.array([1,1,1,1])
y = np.array([0,1,2,3])
f = 1 
ini_params = 1
calc_cross_validation(x,y,f,ini_params)

'''