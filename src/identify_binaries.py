#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
import os
from datetime import datetime

def add_binary_column(catalog):
    log_file = open('log.txt','a')
    log_file.write('\n')
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    log_file.write("Log's date and time: {}\n".format(dt_string))
    log_file.flush()
    
    if(os.path.exists('Catalogs/results_binaries/')):
        print('Catalogs/results_binaries/ exists')
    else:
        os.mkdir('Catalogs/results_binaries/')
    source_id_binaries,total = remove_potential_binaries(catalog)
    
    source_id = catalog['gaia_source_id']
    potential_binary = np.zeros(len(source_id))
    for i,x in enumerate(source_id):
        if(x in source_id_binaries):
            potential_binary[i] = 1
    catalog['potential_binary'] = potential_binary
    n_b_tot = len(source_id_binaries)
    n_b = len(set(source_id_binaries))
    log_file.write('Number of photometric binaries: {}\n'.format(n_b_tot))
    log_file.write('Number of single stars: {}\n'.format(n_b))
    return catalog
    
def get_binaries(label,model,g_rp,g_abs,source_id):
    N = 30
    g_rp_limit = 1.2
    mask_binaries = ((abs(np.polyval(np.flipud(model['cf']),g_rp) 
                          - (g_abs+0.75)) < 0.2) 
                     * (g_rp<g_rp_limit))

    x_group = np.linspace(model['color_range'][0],model['color_range'][1],N)
    y_group = np.polyval(np.flipud(model['cf']),x_group)

    plt.scatter(g_rp,g_abs,label=label)
    plt.scatter(g_rp[mask_binaries],g_abs[mask_binaries],facecolor='none',
                edgecolor='k',s=80,label='binaries')
    plt.plot(x_group,y_group,'--k')
    plt.plot(x_group,y_group-0.75,'--r')
    plt.xlabel('g-rp')
    plt.ylabel('g_abs')
    plt.legend()
    plt.ylim(14,4)
    plt.xlim(0.6,1.5)
    plt.savefig('Catalogs/results_binaries/cmd_'+label)
    plt.close()

    return source_id[mask_binaries]

def remove_potential_binaries(catalog):
    parallax = catalog['parallax']
    g = catalog['phot_g_mean_mag']
    rp = catalog['phot_rp_mean_mag']
    g_rp = g-rp
    g_abs = catalog['phot_g_mean_mag'] - 5*(np.log10(10**3/parallax)-1)
    group_name = catalog['group_name']
    source_id = catalog['gaia_source_id']
    
    path = 'Catalogs/Models/moving_groups/'
    PLE_model = readsav(path+'dr2_seq_fit_cf_PLE.sav')
    HYA_model = readsav(path+'dr2_seq_fit_cf_HYA.sav')
    CBER_model = readsav(path+'dr2_seq_fit_cf_CBER.sav')
    
    mask_abdmg = (group_name=='ABDMG')
    mask_cber = (group_name=='CBER') 
    mask_hya = group_name=='HYA'
    
    models = [PLE_model,HYA_model,CBER_model]
    masks_groups = [mask_abdmg,mask_hya,mask_cber]
    labels = ['PLE','HYA','CBER']
    
    source_id_binaries = []
    total = 0
    
    for model_i,mask_group_i,label_i in zip(models,masks_groups,labels):
        gaia_id = get_binaries(label_i,model_i,g_rp[mask_group_i],
                               g_abs[mask_group_i],source_id[mask_group_i])
        for x in gaia_id:
            source_id_binaries.append(x)
        total+=len(g_rp[mask_group_i])
    
    return source_id_binaries,total