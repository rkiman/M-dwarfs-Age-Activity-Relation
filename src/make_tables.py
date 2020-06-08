#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from astropy.table import Table
import numpy as np
from .accretors import calc_delta_ha_for_accretors
    

def make_summary_sources(source_num,Ncomp,compatible,total_comp,order,
                                 overlap_not_comp,total_overlap_comp,
                                 total_overlap_not_comp):

    #Make latex table summarizing all they studies used in this paper. 

     
    #Set up path to save file
    dat_path='data/'
    if(os.path.exists(dat_path + 'sources_summary.txt')):
        os.remove(dat_path + 'sources_summary.txt')
    file_sources = open(dat_path + 'sources_summary.txt','x')
    file_sources.write('#Ref Nstars Noverlap Ncompatible OC\n')
                       
    #Load reference table for source
    source_ref_table = Table.read('data/source_ref.csv')
    n_sources = len(source_ref_table['source_ref'])
    
    #Identify number for Kiman et al. 2019
    mask_kiman = source_ref_table['source_ref']=='Kiman 2019'
    i_kiman = source_ref_table['source_num'][mask_kiman][0]


    for i in range(n_sources):  
        #Total number of objects from the source
        n_i = len(source_num[source_num == i]) 
        
        #Sources with objects compatible objects from the soruce
        if i in compatible:
            mask = i == compatible
            n_compatible_i = int(total_comp[mask][0])
            n_overlap_i = int(total_overlap_comp[mask][0])
            if i==i_kiman:
                order_i = 0
            else:
                order_i = int(order[mask][0])
            
            file_sources.write(str(source_ref_table['source_num'][i]) + '\t' + 
                               str(n_i) +'\t'+ str(n_overlap_i) + '\t' + 
                               str(n_compatible_i) + '\t' + 
                               str(order_i) + '\n') 
        #Sources with objects not compatible but with overlap
        elif(i in overlap_not_comp):
            mask = i == overlap_not_comp
            n_compatible_i = 0
            order_i = np.nan
            n_overlap_i = int(total_overlap_not_comp[mask][0])
            file_sources.write(str(source_ref_table['source_num'][i]) + '\t' + 
                           str(n_i) +'\t'+ str(n_overlap_i) + '\t' + 
                           str(n_compatible_i) + '\t' + 
                           str(order_i) + '\n') 
        #Sources with objects without overlap
        else:
            n_compatible_i = 0
            order_i = np.nan
            n_overlap_i = 0
            file_sources.write(str(source_ref_table['source_num'][i]) + '\t' + 
                           str(n_i) +'\t'+ str(n_overlap_i) + '\t' + 
                           str(n_compatible_i) + '\t' + 
                           str(order_i) + '\n') 

    return 0


def make_table_for_wd(binaries):
    '''
    Makes latex table with the summary of the white dwarfs values calculated
    with wdwarfdate.
    '''
    
    dropbox_path = '/Users/rociokiman/Dropbox/Apps/Overleaf/'
    paper = 'wdwarfdate/'
    os.remove(dropbox_path + paper + 'wd_summary.tex')
    file_sources = open(dropbox_path + paper + 'wd_summary.tex','x')
    
    n_wd = 5
    rand_idx = np.random.randint(0,len(binaries),n_wd)

    #Header
    file_sources.write('\\begin{deluxetable*}{ccccccc}[ht!]\n')
    file_sources.write('\\tablewidth{290pt}\n')
    file_sources.write('\\tabletypesize{\scriptsize}\n')
    file_sources.write('\\tablecaption{WD sample. \\label{table:wd_sample}}\n')
    file_sources.write('\\tablehead{\
                        \\colhead{$T_{\\rm eff}$}\
                        & \\colhead{$\log (g)$}\
                        & \\colhead{$t_{\\rm ms}$(Gyr)}\
                        & \\colhead{$t_{\\rm cool}$(Gyr)}\
                        & \\colhead{$t_{\\rm tot}$(Gyr)}\
                        & \\colhead{$m_{\\rm i}$(M\\textsubscript{\(\odot\)})}\
                        & \\colhead{$m_{\\rm f}$(M\\textsubscript{\(\odot\)})}\
                        \n}')
    file_sources.write('\\startdata \n')
    
    for i in rand_idx:
        #teff
        val = str(np.round(binaries['TeffH'][i],2))
        err = str(np.round(binaries['e_TeffH'][i],2))
        teff = '$' + val + '\pm ' + err + '$'
        #logg
        val = str(np.round(binaries['loggH'][i],2))
        err = str(np.round(binaries['e_loggH'][i],2))
        logg = '$' + val + '\pm ' + err + '$'
        #ms age
        val = str(np.round(binaries['ms_age_median'][i]/1e9,2))
        err_low = str(np.round(binaries['ms_age_err_low'][i]/1e9,2))
        err_high = str(np.round(binaries['ms_age_err_high'][i]/1e9,2))
        ms = '$' + val + '_{-' + err_low + '}^{+' + err_high + '}$'
        #cool age
        val = str(np.round(binaries['cooling_age_median'][i]/1e9,2))
        err_low = str(np.round(binaries['cooling_age_err_low'][i]/1e9,2))
        err_high = str(np.round(binaries['cooling_age_err_high'][i]/1e9,2))
        cool = '$' + val + '_{-' + err_low + '}^{+' + err_high + '}$'
        #tot age
        val = str(np.round(binaries['total_age_median'][i]/1e9,2))
        err_low = str(np.round(binaries['total_age_err_low'][i]/1e9,2))
        err_high = str(np.round(binaries['total_age_err_high'][i]/1e9,2))
        tot = '$' + val + '_{-' + err_low + '}^{+' + err_high + '}$'
        #initial mass
        val = str(np.round(binaries['initial_mass_median'][i],2))
        err_low = str(np.round(binaries['initial_mass_err_low'][i],2))
        err_high = str(np.round(binaries['initial_mass_err_high'][i],2))
        mini = '$' + val + '_{-' + err_low + '}^{+' + err_high + '}$'
        #final mass
        val = str(np.round(binaries['final_mass_median'][i],2))
        err_low = str(np.round(binaries['final_mass_err_low'][i],2))
        err_high = str(np.round(binaries['final_mass_err_high'][i],2))
        mfin = '$' + val + '_{-' + err_low + '}^{+' + err_high + '}$'
        #complete line
        file_sources.write(teff+'&'+logg+'&'+ms+'&'+cool+'&'+tot+'&'
                           +mini+'&'+mfin+'\\\ \n') 

    file_sources.write('\\enddata \n')
    file_sources.write('\\end{deluxetable*}\n')

    return 0

def make_table_wd_ages(binaries):
    '''
    Makes a table summary of white dwarfs with only the total age calculated
    and the Gaia id for the M dwarf and the white dwarf in the pair.
    '''

    dropbox_path = '/Users/rociokiman/Dropbox/Apps/Overleaf/'
    paper = 'Age-Activity Relation for M dwarfs/'
    if(os.path.exists(dropbox_path + paper + 'wd_ages.tex')):
        os.remove(dropbox_path + paper + 'wd_ages.tex')
    file_sources = open(dropbox_path + paper + 'wd_ages.tex','x')
    
    n_wd = 5
    rand_idx = np.random.randint(0,len(binaries),n_wd)

    #Header
    file_sources.write('\\begin{deluxetable*}{ccc}[ht!]\n')
    file_sources.write('\\tablewidth{290pt}\n')
    file_sources.write('\\tabletypesize{\scriptsize}\n')
    file_sources.write('\\tablecaption{WD sample. \\label{table:wd_sample}}\n')
    file_sources.write('\\tablehead{\
                        \\multicolumn{2}{c}{\\textit{Gaia} source id} \
                        & \\colhead{Total age} \
                        \\\ M dwarf & White dwarf & (Gyr) \n}')
    file_sources.write('\\startdata \n')
    
    for i in rand_idx:
        source_id_m = '$' + str(np.round(binaries['Source'][i],2)) + '$'
        source_id_wd = '$' + str(np.round(binaries['m_source_id'][i],2)) + '$'
        age = str(np.round(binaries['total_age_median'][i]/1e9,2))
        err_low = str(np.round(binaries['total_age_err_low'][i]/1e9,2))
        err_high = str(np.round(binaries['total_age_err_high'][i]/1e9,2))
        tot = '$' + age + '_{-' + err_low + '}^{+' + err_high + '}$'
        file_sources.write(source_id_m+'&'+source_id_wd+'&'+tot+'\\\ \n') 

    file_sources.write('\\enddata \n')
    file_sources.write('\\end{deluxetable*}\n')

    return 0

def make_table_summary_age_calibrators(age_calibrators):
    '''
    Makes a table summarizing the age calibrators and the sources from where
    they came from.
    '''

    mask_ha = ~np.isnan(age_calibrators['ewha'])
    source_ref_table = Table.read('data/source_ref.csv')
    data_compatible = np.loadtxt('data/sources_summary.txt')

    dropbox_path = '/Users/rociokiman/Dropbox/Apps/Overleaf/'
    paper = 'Age-Activity Relation for M dwarfs/'
    if(os.path.exists(dropbox_path + paper + 'summary_age_calibrators.tex')):
        os.remove(dropbox_path + paper + 'summary_age_calibrators.tex')
    file_sources = open(dropbox_path + paper + 'summary_age_calibrators.tex',
                        'x')
    
    #Header
    header_text = '\\tablehead{\
    \\colhead{Reference \\tablenotemark{a}} \
    & \\colhead{Spectral} \
    & \\multicolumn{2}{c}{$N$ of M dwarfs} \
    & \\colhead{OC \\tablenotemark{c}} \
    & \\multicolumn{2}{c}{Ages from} \
    \\\ & Resolution & Total & Compatible \\tablenotemark{b} \
    &  &moving group & white dwarf \n}'
    
    title = '\\tablecaption{Age Calibrators summary. \\label{table:age_cal}}\n'
    
    file_sources.write('\\begin{deluxetable*}{lcccccc}[ht!]\n')
    file_sources.write('\\tablewidth{290pt}\n')
    file_sources.write('\\tabletypesize{\scriptsize}\n')
    file_sources.write(title)
    file_sources.write(header_text)
    file_sources.write('\\startdata \n')
    
    catalogs_no_overlap = []
    catalogs_no_compatible = []
    other_catalogs = []
    catalogs_no_age = []
    
    for x,y,z in zip(source_ref_table['source_num'],source_ref_table['cite'],
                     source_ref_table['resolution']):
        mask_source = age_calibrators['source_num'] == x
        mask_source_2 = data_compatible[:,0] == x
        mask_wd = age_calibrators['group_num'] == 0
        mask_mg = age_calibrators['group_num'] != 0
        n_wd = len(age_calibrators[mask_source*mask_wd*mask_ha])
        n_mg = len(age_calibrators[mask_source*mask_mg*mask_ha])
        n_wd_all = len(age_calibrators[mask_source*mask_wd])
        n_mg_all = len(age_calibrators[mask_source*mask_mg])
        
        n_tot = int(data_compatible[:,1][mask_source_2])
        n_overlap = int(data_compatible[:,2][mask_source_2])
        n_comp = int(data_compatible[:,3][mask_source_2])
        if(~np.isnan(data_compatible[:,4][mask_source_2])):
            n_oc = int(data_compatible[:,4][mask_source_2])
        else:
            n_oc = data_compatible[:,4][mask_source_2]
        
        if(n_tot==0 and n_overlap==0 and n_comp==0):
            other_catalogs.append(y)
        elif(n_tot!=0 and n_overlap==0 and n_comp==0):
            catalogs_no_overlap.append(y)
        elif(n_tot!=0 and n_overlap!=0 and n_comp==0):
            catalogs_no_compatible.append(y)
        else:
            if(any(np.array([n_wd,n_mg,n_wd_all,n_mg_all])!=0)):
                if(all(np.array([n_wd,n_wd_all])==0)):
                    if(n_mg==n_mg_all):
                        file_sources.write(y+'&'+z+'&'+str(n_tot)+'&'
                                           +str(n_comp)+
                                           '&'+str(n_oc)+'&'+str(n_mg)+
                                           '&'+'-'+'\\\ \n')
                    else:
                        file_sources.write(y+'&'+z+'&'+str(n_tot)+'&'
                                           +str(n_comp)+'&'
                                           +str(n_oc)+'&'+str(n_mg)+'/'
                                           +str(n_mg_all)+'&'
                                           +'-'+'\\\ \n')
                elif(all(np.array([n_mg,n_mg_all])==0)):
                    if(n_wd==n_wd_all):
                        file_sources.write(y+'&'+z+'&'+str(n_tot)+'&'
                                           +str(n_comp)+
                                           '&'+str(n_oc)+'&'+'-'+
                                           '&'+str(n_wd)+'\\\ \n')
                    else:
                        file_sources.write(y+'&'+z+'&'+str(n_tot)+'&'
                                           +str(n_comp)+'&'
                                           +str(n_oc)+'&'+'-'+'&'
                                           +str(n_wd)+'/'+str(n_wd_all)
                                           +'\\\ \n')
                elif(n_mg==n_mg_all and n_wd==n_wd_all):
                    file_sources.write(y+'&'+z+'&'+str(n_tot)+'&'
                                       +str(n_comp)+'&'
                                       +str(n_oc)+'&'+str(n_mg)+'&'
                                       +str(n_wd)+'\\\ \n') 
                else:
                    file_sources.write(y+'&'+z+'&'+str(n_tot)+'&'+str(n_comp)
                    +'&'+str(n_oc)+'&'+str(n_mg)+'/'+str(n_mg_all)+'&'
                                       +str(n_wd)+'/'+str(n_wd_all)+'\\\ \n')
            else:
                catalogs_no_age.append(y)

    file_sources.write('\\enddata \n')

    papers_notes='Compatible catalogs without age calibrators: '
    for catalog_i in catalogs_no_age:
        papers_notes = papers_notes + catalog_i + ', '
        
    papers_notes=papers_notes+'.\n Catalogs with overlap but not compatibles: '
    for catalog_i in catalogs_no_compatible:
        papers_notes = papers_notes + catalog_i + ', '

    papers_notes=papers_notes+'.\n Catalogs without overlap: '
    for catalog_i in catalogs_no_overlap:
        papers_notes = papers_notes + catalog_i + ', '
    
    papers_notes=papers_notes+'.\n Other catalogs checked: '
    for catalog_i in other_catalogs:
        papers_notes = papers_notes + catalog_i + ', '
    
    file_sources.write('\\tablenotetext{a}{'+ papers_notes +'.}\n')
    file_sources.write('\\tablenotetext{b}{Compatible with \
                       \\citet{Kiman2019}.}\n')
    file_sources.write('\\tablenotetext{c}{Order of compatibility. Order $1$\
                       is compatible with \\citet{Kiman2019}. Order $2$ is \
                       compatible with at least one order $1$ catalog.}\n')
    file_sources.write('\\end{deluxetable*}\n')

    return 0

def make_table_accretors(accretors):
    '''
    Makes a latex table with a summary of the accreators found in the
    literature search sample using the criterion from 
    White,R.J. & Basri,G. Astrophys. J. 582, 1109–1122 (2003).
    '''
    mask_nan = ~np.isnan(accretors['ewha_error'])
    n_acc = 5
    delta_ha = np.ones(n_acc)*-1
    while(any(delta_ha<0)):
        rand_idx = np.random.randint(0,len(accretors[mask_nan]),n_acc)    
        source_id = accretors['source_id'][mask_nan][rand_idx]
        g = accretors['g_corr'][mask_nan][rand_idx]
        rp = accretors['rp_corr'][mask_nan][rand_idx]
        g_rp = g-rp
        ewha = accretors['ewha'][mask_nan][rand_idx]
        ewha_error = accretors['ewha_error'][mask_nan][rand_idx]
        d = [calc_delta_ha_for_accretors(x,y) for x,y in zip(g_rp,ewha)]
        delta_ha = np.array([np.round(x,2) for x in d])
    
    
    dropbox_path = '/Users/rociokiman/Dropbox/Apps/Overleaf/'
    paper = 'Age-Activity Relation for M dwarfs/'
    if(os.path.exists(dropbox_path + paper + 'summary_accretors.tex')):
        os.remove(dropbox_path + paper + 'summary_accretors.tex')
    file_sources = open(dropbox_path + paper + 'summary_accretors.tex','x')
    
    #Header
    header_text = '\\tablehead{\
    \\colhead{\\textit{Gaia} source id} \
    & \\colhead{$\\haew$} \
    & \\colhead{$\\Delta \\haew$\\tablenotemark{a}} \n}'
    
    title = '\\tablecaption{Accretors found with the criterion from \
    \\citet{White2003}. \\label{table:accretors}}\n'
    
    file_sources.write('\\begin{deluxetable*}{lcc}[ht!]\n')
    file_sources.write('\\tablewidth{290pt}\n')
    file_sources.write('\\tabletypesize{\scriptsize}\n')
    file_sources.write(title)
    file_sources.write(header_text)
    file_sources.write('\\startdata \n')
    
    
    for x,y,z,w in zip(source_id,ewha,ewha_error,delta_ha):
        file_sources.write(str(x)+'& $'+ 
                           str(np.round(y,2)) + '\\pm' 
                           +str(np.round(z,2)) +'$ & $'
                           +str(w) +'$ \\\ \n')
    file_sources.write('\\enddata \n')
    file_sources.write('\\tablenotetext{a}{Delta above the $\\haew$ limit.}\n')
    file_sources.write('\\end{deluxetable*}\n')
    
