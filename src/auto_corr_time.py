#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
The code for this file comes from the tutorial by Dan Foreman-Mackey:
https://emcee.readthedocs.io/en/stable/tutorials/autocorr/
'''

import numpy as np

def calc_auto_corr_time(chain):

    N = np.exp(np.linspace(np.log(100),np.log(chain.shape[1]), 10)).astype(int)
    corr_time = []
    #Loop over parameters
    for j in range(0,3):
        chain1 = chain[:, :, j]
        new = np.empty(len(N))
        for i, n in enumerate(N):
            auto_corr = autocorr_new(chain1[:, :n])
            new[i] = auto_corr
        corr_time.append(new[-1])

    return np.array(corr_time)

def next_pow_two_wd(n):
    i = 1
    while i < n:
        i = i << 1
    return i

def autocorr_func_1d_wd(x, norm=True):
    x = np.atleast_1d(x)
    if len(x.shape) != 1:
        raise ValueError("invalid dimensions for 1D autocorrelation function")
    n = next_pow_two_wd(len(x))

    # Compute the FFT and then (from that) the auto-correlation function
    f = np.fft.fft(x - np.mean(x), n=2 * n)
    acf = np.fft.ifft(f * np.conjugate(f))[: len(x)].real
    acf /= 4 * n

    # Optionally normalize
    if norm:
        acf /= acf[0]

    return acf

def auto_window_wd(taus, c):
    m = np.arange(len(taus)) < c * taus
    if np.any(m):
        return np.argmin(m)
    return len(taus) - 1

def autocorr_new(y, c=5.0):
    f = np.zeros(y.shape[1])
    for yy in y:
        f += autocorr_func_1d_wd(yy)
    f /= len(y)
    taus = 2.0 * np.cumsum(f) - 1.0
    window = auto_window_wd(taus, c)
    return taus[window]