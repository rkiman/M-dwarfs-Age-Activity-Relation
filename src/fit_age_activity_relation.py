#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy.optimize as op
import emcee
import corner
import matplotlib.pyplot as plt
from .auto_corr_time import calc_autocorrelation_time

def west2008(x,a,b,n,l):
    lhalbol = np.ones(len(x))*np.nan
    mask = x < l
    lhalbol[mask] = a/(x[mask]**n-l**n) - b
    return lhalbol

def fit_halpha_bpl(params,log_age):

    a0,a1,a2,a3 = params
    
    halpha_model = np.ones(len(log_age))*np.nan
    mask = log_age < a0
    halpha_model[mask] = a1*log_age[mask] + a2
    halpha_model[~mask] = a3*log_age[~mask] + (a1-a3)*a0 + a2 

    return halpha_model
        
def lnlike_age_bpl(params,log_age,log_lhalbol,log_lhalbol_error):
    
    a0,a1,a2,a3 = params
    
    model_halpha = fit_halpha_bpl(params,log_age)
    sigma2 = log_lhalbol_error ** 2 
    if(a0<0 or a0>10.3 or a3>100 or a3<-100 or a1>100 or a1 < -100):
        return -np.inf
    else:
        return -0.5 * np.sum((log_lhalbol - model_halpha) ** 2 / sigma2)



def fit_relation_bpl(mask,log_age,log_lhalbol,
                     log_lhalbol_error,ini_params = np.array([9,-.1,1,-4]),
                     sigma_random = 0.001,
                     name='corner_fit.png',
                     name_file_auto_corr='auto_corr.png'):

    #Optimize the parameters with scipy. Gets close to the result
    nll = lambda *args: -lnlike_age_bpl(*args)
    params = op.minimize(nll, ini_params, 
                         args=(log_age[mask],log_lhalbol[mask],
                               log_lhalbol_error[mask]))
    
    #New initial parameters are the results from the optimization
    ini_params = params.x
    n_params = len(params.x)
    
    #Set dimension and walkder
    ndim, nwalkers = n_params, 200

    #Initialize every walker
    p0 = np.array([ini_params+(np.random.rand(n_params)-0.5)*sigma_random
                   for i in range(nwalkers)])

    #Define sampler
    sampler = emcee.EnsembleSampler(nwalkers,ndim,lnlike_age_bpl, 
                                    args=[log_age[mask],log_lhalbol[mask],
                                          log_lhalbol_error[mask]])
    
    print('Running burn in')
    p,_,_ = sampler.run_mcmc(p0, 2000)
    
    #Take medians of the positions of each walkers to initialize again
    chain = sampler.chain[:,:,:]
    flat_samples = chain.reshape((-1,ndim))
    p0_medians = np.array([np.median(flat_samples[:,x]) for x in range(n_params)])
    p0_new = np.array([p0_medians+(np.random.rand(n_params)-0.5)*sigma_random
                       for i in range(nwalkers)])
    
    sampler.reset()
    
    print('Running mcmc to calculate auto-correlation time')
    p,_,_ = sampler.run_mcmc(p0_new, 100000)
    chain = sampler.chain
    autoc_time = calc_autocorrelation_time(chain,name_file_auto_corr)
    sampler.reset()
    
    print('Running mcmc to calculate parameters')
    p,_,_ = sampler.run_mcmc(p, 1000*autoc_time)
    chain = sampler.chain[:,:,:]
    flat_samples = chain.reshape((-1,ndim))
    
    plot_corner_plot(flat_samples,name)
    
    return chain,flat_samples

def plot_corner_plot(flat_samples,name):
    #labels = ['$a_0$','$a_1$','$a_2$','$a_3$','logf']
    labels = ['$a_0$','$a_1$','$a_2$','$a_3$']

    fig = corner.corner(flat_samples,labels=labels,quantiles=[.16,.50,.84],
                        show_titles=True, title_kwargs={"fontsize": 12})
    dropbox = '/Users/rociokiman/Dropbox (Personal)/Apps/Overleaf'
    path = dropbox+'/Age-Activity Relation for M dwarfs/'+name
    fig.savefig(path,dpi=300)
    return 0
    
    
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


'''
def fit_halpha(params,color,log_age):

    a0,a1,a2,a3,a4,a5 = params
    
    halpha_model = (a0 + a1*color + a2*color**2 +
                    a3*log_age + a4*log_age**2 +
                    a5*log_age*color)

    return halpha_model
        
def lnlike_color_age(params,color,log_age,log_lhalbol,log_lhalbol_error):
    model_halpha = fit_halpha(params,color,log_age)
    return -0.5*np.sum((model_halpha - log_lhalbol)**2/(log_lhalbol_error**2)) 


def fit_relation_complex_func(mask,color,log_age,log_lhalbol,
                              log_lhalbol_error,name='corner_fit.png'):

    ini_params = np.ones(6)*0.1
    
    nll = lambda *args: -lnlike_color_age(*args)
    params = op.minimize(nll, ini_params, 
                         args=(color[mask],log_age[mask],log_lhalbol[mask],
                               log_lhalbol_error[mask]))
    
    n_params = len(params.x)
    ndim, nwalkers = n_params, 100
    
    p0 = np.array([params.x+np.random.rand(n_params)*2-1
                   for i in range(nwalkers)])
    
    sampler = emcee.EnsembleSampler(nwalkers,ndim,lnlike_color_age, 
                                    args=[color[mask],log_age[mask],
                                          log_lhalbol[mask],
                                          log_lhalbol_error[mask]])
    sampler.run_mcmc(p0, 10000)
    
    chain = sampler.chain[:,500:,:]
    flat_samples = chain.reshape((-1,ndim))
    
    labels = ['$a_0$','$a_1$','$a_2$','$a_3$','$a_4$','$a_5$']

    fig = corner.corner(flat_samples,labels=labels,quantiles=[.16,.50,.84],
                        show_titles=True, title_kwargs={"fontsize": 12})
    dropbox = '/Users/rociokiman/Dropbox (Personal)/Apps/Overleaf'
    path = dropbox+'/Age-Activity Relation for M dwarfs/'+name
    fig.savefig(path,dpi=300)
    
    return flat_samples

    

def lnlike_age_bpl_var(params,log_age,log_lhalbol,log_lhalbol_error):
    
    a0,a1,a2,a3,log_f = params
    
    model_halpha = fit_halpha_bpl(np.array([a0,a1,a2,a3]),log_age)
    sigma2 = log_lhalbol_error ** 2 + model_halpha ** 2 * np.exp(2 * log_f)
    
    if(a0<0 or a0>10.3 or a3>0 or a3<-100 or a1>100 or a1 < -100
       or a2 < 0 or -10.0 > log_f or log_f > 1.0):
        return -np.inf
    else:
        return -0.5 * np.sum((log_lhalbol - model_halpha) ** 2 / sigma2 + np.log(sigma2))

    
def fit_relation_bpl_var(mask,log_age,log_lhalbol,
                         log_lhalbol_error,name='corner_fit.png'):

    #ini_params = np.array([9,-.1,1,-4,.01])
    ini_params = np.array([9,-.1,1,-4,-3])
    
    nll = lambda *args: -lnlike_age_bpl_var(*args)
    params = op.minimize(nll, ini_params, 
                         args=(log_age[mask],log_lhalbol[mask],
                               log_lhalbol_error[mask]))
    
    #ini_params = [9,-1,-4,1,1]
    ini_params = params.x
    print(ini_params)
    n_params = len(params.x)
    
    ndim, nwalkers = n_params, 150

    p0 = np.array([ini_params+np.random.rand(n_params)*0.0001
                   for i in range(nwalkers)])
    
    sampler = emcee.EnsembleSampler(nwalkers,ndim,lnlike_age_bpl_var, 
                                    args=[log_age[mask],log_lhalbol[mask],
                                          log_lhalbol_error[mask]])
    sampler.run_mcmc(p0, 10000)
    
    chain = sampler.chain[:,500:,:]
    flat_samples = chain.reshape((-1,ndim))
    
    #labels = ['$a_0$','$a_1$','$a_2$','$a_3$','logf']
    labels = ['$a_0$','$a_1$','$a_2$','$a_3$',r'$\log f$']

    fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
    for i in range(ndim):
        ax = axes[i]
        ax.plot(flat_samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(flat_samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)
    
    axes[-1].set_xlabel("step number");

    fig = corner.corner(flat_samples,labels=labels,quantiles=[.16,.50,.84],
                        show_titles=True, title_kwargs={"fontsize": 12})
    dropbox = '/Users/rociokiman/Dropbox (Personal)/Apps/Overleaf'
    path = dropbox+'/Age-Activity Relation for M dwarfs/'+name
    fig.savefig(path,dpi=300)
    
    return chain,flat_samples
'''
#def plot_walkers()
