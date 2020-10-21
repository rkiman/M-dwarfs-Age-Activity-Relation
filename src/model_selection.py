#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from .cross_validation import calc_cross_validation
from scipy.optimize import curve_fit

def linear_model(x,a,b):
    return a*x+b

def polynomial(x,*p):
    return np.polyval(p,x)

def broken_power_law(x,*p):
    log_t0,alpha1,beta1,alpha2,_ = p
    if(isinstance(x, float)):
        if(x<log_t0):
            y_model = alpha1*(x-log_t0) + beta1
        else:
            y_model = alpha2*(x-log_t0) + beta1
    else:
        y_model = np.ones(len(x))*np.nan
        mask = x < log_t0
        y_model[mask] = alpha1*(x[mask]-log_t0) + beta1
        y_model[~mask] = alpha2*(x[~mask]-log_t0) + beta1
    
    return y_model


def select_model(lage,log_lhalbol,log_lhalbol_error,mask_good,
                 bpl_ini_p=np.array([9,-.01,-3,-4])):
    
    scores = []
    x_degree = np.arange(2,14)
    for n in x_degree:
        ini_params = np.ones(n)
        scores.append(calc_cross_validation(lage[mask_good],
                                            log_lhalbol[mask_good],
                                            log_lhalbol_error[mask_good],
                                            polynomial,ini_params))
        
    score_broken = calc_cross_validation(lage[mask_good],
                                         log_lhalbol[mask_good],
                                         log_lhalbol_error[mask_good],
                                         broken_power_law,bpl_ini_p)

    
    ini_params = np.ones(6)
    popt,pcov=curve_fit(polynomial,lage[mask_good],log_lhalbol[mask_good],
                        p0=ini_params,
                        sigma=log_lhalbol_error[mask_good])
    
    ini_params = np.ones(2)
    popt_linear,pcov=curve_fit(polynomial,lage[mask_good],
                               log_lhalbol[mask_good],
                               p0=ini_params,
                               sigma=log_lhalbol_error[mask_good])
    
    ini_params = bpl_ini_p
    popt_broken,pcov=curve_fit(broken_power_law,lage[mask_good],
                               log_lhalbol[mask_good],
                               p0=ini_params,
                               sigma=log_lhalbol_error[mask_good])
    
    x=np.linspace(6,10,100)
    f,(ax1,ax2) = plt.subplots(1,2,figsize=(7,3))
    ax1.plot(x_degree-1,scores,'.-',label='Polynomial')
    ax1.axhline(y=score_broken,color='k',label='Broken power law')
    ax1.set_ylabel('Cross-validation score')
    ax1.set_xlabel('Polynomial degree')
    ax1.legend()

    ax2.errorbar(lage[mask_good],log_lhalbol[mask_good],
                 yerr=log_lhalbol_error[mask_good],
                 fmt='.',elinewidth=.5)
    ax2.plot(x,np.polyval(popt,x),color='r',label='Polynomial of 3rd degree')
    ax2.plot(x,np.polyval(popt_linear,x),color='purple',label='Linear fit')
    ax2.plot(x,broken_power_law(x,*popt_broken),color='g',
             label='Broken power law')
    ax2.set_xlabel(r'$\log _{10}Age/$yr')
    ax2.set_ylabel(r'$\log _{10}L_{\rm H\alpha}/L_{\rm BOL}$')
    ax2.legend()
    plt.tight_layout()
    plt.show() 
    
    return scores,score_broken
    
