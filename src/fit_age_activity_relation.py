#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import emcee
import matplotlib.pyplot as plt
from .auto_corr_time import calc_auto_corr_time

def west2008(x,a,b,n,l):
    lhalbol = np.ones(len(x))*np.nan
    mask = x < l
    lhalbol[mask] = a/(x[mask]**n-l**n) - b
    return lhalbol

def fit_halpha_bpl(params,log_age):

    log_t0,alpha1,beta1,alpha2,_ = params
    
    halpha_model = np.ones(len(log_age))*np.nan
    mask = log_age < log_t0
    halpha_model[mask] = alpha1*(log_age[mask]-log_t0) + beta1
    halpha_model[~mask] = alpha2*(log_age[~mask]-log_t0) + beta1

    return halpha_model

    
def lnlike_age_bpl(params,log_age,log_lhalbol,log_lhalbol_error):
    
    log_t0,alpha1,beta1,alpha2,sigma_v = params
    
    model_halpha = fit_halpha_bpl(params,log_age)
    sigmalpha2 = log_lhalbol_error ** 2 + sigma_v**2
    if(log_t0<0 or log_t0>10.3 or alpha2>100 or alpha2<-100 or 
       alpha1>100 or alpha1 < -100 
       or sigma_v < 0 or sigma_v > 10):
        return -np.inf
    else:
        return -0.5 * np.sum((log_lhalbol - model_halpha) ** 2 / sigmalpha2 
                             + np.log(sigmalpha2))



def fit_relation_bpl(mask,log_age,log_lhalbol,
                     log_lhalbol_error,ini_params = np.array([9,-.1,1,-4,1]),
                     sigma_random = 0.001,
                     name_file_auto_corr='auto_corr.png'):
    max_n = 100000
    n_indep_samples = 100
    #Optimize the parameters with scipy. Gets close to the result
    #nll = lambda *args: -lnlike_age_bpl(*args)
    #params = op.minimize(nll, ini_params, 
    #                     args=(log_age[mask],log_lhalbol[mask],
    #                           log_lhalbol_error[mask]))
    
    #New initial parameters are the results from the optimization
    #ini_params = params.x
    #print(ini_params)
    n_params = len(ini_params)
    
    #Set dimension and walkder
    ndim, nwalkers = n_params, 200

    #Initialize every walker
    p0 = np.array([ini_params+(np.random.rand(n_params)-0.5)*sigma_random
                   for i in range(nwalkers)])

    #Define sampler
    sampler = emcee.EnsembleSampler(nwalkers,ndim,lnlike_age_bpl, 
                                    args=[log_age[mask],log_lhalbol[mask],
                                          log_lhalbol_error[mask]])
    #run burn in
    p0_new,_,_ = sampler.run_mcmc(p0, 1000)
    sampler.reset()
    
    n_steps = int(max_n/100)
    
    index = 0
    autocorr = np.empty(max_n)
    
    # This will be useful to testing convergence
    old_tau = np.inf
    
    # going to run the mcmc in groups of 100 steps
    #The code for the calculation of the autocorrelation time is from the
    #tutorial in https://emcee.readthedocs.io/en/stable/tutorials/autocorr/
    for x in range(n_steps):
        p0_new,_,_ = sampler.run_mcmc(p0_new, 100)
        chain = sampler.chain
        # Compute the autocorrelation time so far
        tau = calc_auto_corr_time(chain)
        autocorr[index] = np.mean(tau)
        index += 1
    
        # Check convergence
        converged = np.all(tau * n_indep_samples < (x+1)*100)
        converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
        if converged:
            break
        old_tau = tau
        
    N = 100 * np.arange(1, index + 1)
    plt.plot(N, N / 100.0, "--k", label=r"$\tau = N/100$")
    plt.loglog(N, autocorr[:index], "-")
    plt.xlabel("number of samples, $N$")
    plt.ylabel(r"mean $\hat{\tau}$")
    plt.legend(fontsize=14)
    plt.grid()
    plt.savefig(name_file_auto_corr)
    plt.close()
    
    chain = sampler.chain[:,:,:]
    flat_samples = chain.reshape((-1,ndim))
    
    return chain,flat_samples
    
    
def plot_walkers(chain_simple,labels=[]):
    n_params = len(chain_simple[0,0,:])
    walkers = len(chain_simple[:,0,0])
    fig, axs = plt.subplots(n_params,1, figsize=(6, 6), facecolor='w', 
                            edgecolor='k')
    fig.subplots_adjust(hspace = .5, wspace=.005)
    axs = axs.ravel()
    for i in range(n_params):
        for j in range(walkers):
            axs[i].plot(chain_simple[j,:,i],'k-',linewidth=.01)
        if(len(labels)>0):
            axs[i].set_ylabel(labels[i])
    plt.tight_layout()
    plt.show()


