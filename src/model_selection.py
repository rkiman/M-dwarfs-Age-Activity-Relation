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
    a0,a1,a2,a3 = p
    if(isinstance(x, float)):
        if(x<a0):
            y_model = a1*x + a2
        else:
            y_model = a3*x + (a1-a3)*a0 + a2 
    else:
        y_model = np.ones(len(x))*np.nan
        mask = x < a0
        y_model[mask] = a1*x[mask] + a2
        y_model[~mask] = a3*x[~mask] + (a1-a3)*a0 + a2 
    return y_model



def select_model(lage,log_lhalbol,log_lhalbol_error,mask_good):
    
    scores = []
    for n in range(2,8):
        ini_params = np.ones(n)
        scores.append(calc_cross_validation(lage[mask_good],
                                            log_lhalbol[mask_good],
                                            log_lhalbol_error[mask_good],
                                            polynomial,ini_params))
        
    ini_params = np.array([9,-.1,1,-4])
    score_broken = calc_cross_validation(lage[mask_good],
                                         log_lhalbol[mask_good],
                                         log_lhalbol_error[mask_good],
                                         broken_power_law,ini_params)


    ini_params = np.ones(4)
    popt,pcov=curve_fit(polynomial,lage[mask_good],log_lhalbol[mask_good],
                        p0=ini_params,sigma=log_lhalbol_error[mask_good])
    
    ini_params = np.ones(2)
    popt_linear,pcov=curve_fit(polynomial,lage[mask_good],
                               log_lhalbol[mask_good],
                               p0=ini_params,
                               sigma=log_lhalbol_error[mask_good])
    
    ini_params = np.array([9,-.1,1,-4])
    popt_broken,pcov=curve_fit(broken_power_law,lage[mask_good],
                               log_lhalbol[mask_good],
                               p0=ini_params,
                               sigma=log_lhalbol_error[mask_good])
    
    x=np.linspace(6,10,100)
    f,(ax1,ax2) = plt.subplots(1,2,figsize=(11,4.5))
    ax1.plot(np.arange(2,8)-1,scores,'.-',label='Polynomial')
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
    
