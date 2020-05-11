#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from astropy.table import Table
import numpy as np

def make_table_for_paper_sources(source_num,Ncomp,compatible,total_comp,order):
     
    dropbox_path = '/Users/rociokiman/Dropbox/Apps/Overleaf/'
    paper = 'Age-Activity Relation for M dwarfs/'
    os.remove(dropbox_path + paper + 'sources_summary.tex')
    file_sources = open(dropbox_path + paper + 'sources_summary.tex','x')
    
    source_ref_table = Table.read('data/source_ref.csv')
    
    mask_kiman = source_ref_table['source_ref']=='Kiman 2019'
    i_kiman = source_ref_table['source_num'][mask_kiman][0]
    
    #Ntot = len(ra)
    #Ncomp = len(ls_compatible[~np.isnan(ls_compatible['ewha'])])
    n_sources = len(source_ref_table['source_ref'])
    total_incompatible = 0
    papers_notes = 'Other works checked in literature search but were considered incompatible: '
    #Header
    file_sources.write('\\begin{deluxetable*}{lcccc}[ht!]\n')
    file_sources.write('\\tablewidth{290pt}\n')
    file_sources.write('\\tabletypesize{\scriptsize}\n')
    file_sources.write('\\tablecaption{Overlapping references for $\haew$. \\label{table:source_ref}}\n')
    file_sources.write('\\tablehead{\\colhead{Reference \\tablenotemark{a}} & \\colhead{Spectral} & \\colhead{$N$ Stars} & \\colhead{$N$ Stars} & \\colhead{OC \\tablenotemark{c}} \\\ & Resolution & & Compatible \\tablenotemark{b} &\n}')
    file_sources.write('\\startdata \n')
    
    for i in range(n_sources):          
        if i in compatible:
            mask = i == compatible
            compatible_i = int(total_comp[mask][0])
            if i==i_kiman:
                order_i = 0
            else:
                order_i = int(order[mask][0])
        else:
            compatible_i = 0
            order_i = '-'
        n_i = len(source_num[source_num == i])
        if(n_i!=0 and compatible_i!=0):
            file_sources.write('{}&{}&{}&{}&{}\\\ \n'.format(source_ref_table['cite'][i],source_ref_table['resolution'][i],
                                                   n_i,compatible_i,order_i)) 
        elif(compatible_i==0):
            total_incompatible += n_i
            papers_notes = papers_notes + source_ref_table['cite'][i] + ', '
    file_sources.write('\\enddata \n')
    file_sources.write('\\tablenotetext{a}{'+ papers_notes +'.}\n')
    file_sources.write('\\tablenotetext{b}{Compatible with \\citet{Kiman2019}.}\n')
    file_sources.write('\\tablenotetext{c}{Order of compatibility. Order $1$ is compatible with \\citet{Kiman2019}. Order $2$ is compatible with at least one order $1$ catalog.}\n')
    file_sources.write('\\end{deluxetable*}\n')

    return 0

def make_table_for_wd(binaries):
    
    
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
    file_sources.write('\\tablehead{\\colhead{$T_{\\rm eff}$} & \\colhead{$\log (g)$} & \\colhead{$t_{\\rm ms}$(Gyr)} & \\colhead{$t_{\\rm cool}$(Gyr)} & \\colhead{$t_{\\rm tot}$(Gyr)} & \\colhead{$m_{\\rm i}$(M\\textsubscript{\(\odot\)})} & \\colhead{$m_{\\rm f}$(M\\textsubscript{\(\odot\)})}\n}')
    file_sources.write('\\startdata \n')
    
    for i in rand_idx:
        teff = '$' + str(np.round(binaries['TeffH'][i],2)) + '\pm ' + str(np.round(binaries['e_TeffH'][i],2)) + '$'
        logg = '$' + str(np.round(binaries['loggH'][i],2)) + '\pm ' + str(np.round(binaries['e_loggH'][i],2)) + '$'
        ms = '$' + str(np.round(binaries['ms_age_median'][i]/1e9,2)) + '_{-' + str(np.round(binaries['ms_age_err_low'][i]/1e9,2)) + '}^{+' + str(np.round(binaries['ms_age_err_high'][i]/1e9,2)) + '}$'
        cool = '$' + str(np.round(binaries['cooling_age_median'][i]/1e9,2)) + '_{-' + str(np.round(binaries['cooling_age_err_low'][i]/1e9,2)) + '}^{+' + str(np.round(binaries['cooling_age_err_high'][i]/1e9,2)) + '}$'
        tot = '$' + str(np.round(binaries['total_age_median'][i]/1e9,2)) + '_{-' + str(np.round(binaries['total_age_err_low'][i]/1e9,2)) + '}^{+' + str(np.round(binaries['total_age_err_high'][i]/1e9,2)) + '}$'
        mini = '$' + str(np.round(binaries['initial_mass_median'][i],2)) + '_{-' + str(np.round(binaries['initial_mass_err_low'][i],2)) + '}^{+' + str(np.round(binaries['initial_mass_err_high'][i],2)) + '}$'
        mfin = '$' + str(np.round(binaries['final_mass_median'][i],2)) + '_{-' + str(np.round(binaries['final_mass_err_low'][i],2)) + '}^{+' + str(np.round(binaries['final_mass_err_high'][i],2)) + '}$'
        file_sources.write(teff+'&'+logg+'&'+ms+'&'+cool+'&'+tot+'&'+mini+'&'+mfin+'\\\ \n') 

    file_sources.write('\\enddata \n')
    file_sources.write('\\end{deluxetable*}\n')

    return 0

def make_table_wd_ages(binaries):

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
    file_sources.write('\\tablehead{\\multicolumn{2}{c}{\\textit{Gaia} source id} & \\colhead{Total age} \\\ M dwarf & White dwarf & (Gyr) \n}')
    file_sources.write('\\startdata \n')
    
    for i in rand_idx:
        source_id_m = '$' + str(np.round(binaries['Source'][i],2)) + '$'
        source_id_wd = '$' + str(np.round(binaries['m_source_id'][i],2)) + '$'
        tot = '$' + str(np.round(binaries['total_age_median'][i]/1e9,2)) + '_{-' + str(np.round(binaries['total_age_err_low'][i]/1e9,2)) + '}^{+' + str(np.round(binaries['total_age_err_high'][i]/1e9,2)) + '}$'
        file_sources.write(source_id_m+'&'+source_id_wd+'&'+tot+'\\\ \n') 

    file_sources.write('\\enddata \n')
    file_sources.write('\\end{deluxetable*}\n')

    return 0

def make_table_summary_age_calibrators(age_calibrators):

    mask_ha = ~np.isnan(age_calibrators['ewha'])
    source_ref_table = Table.read('../data/source_ref.csv')
        
    dropbox_path = '/Users/rociokiman/Dropbox/Apps/Overleaf/'
    paper = 'Age-Activity Relation for M dwarfs/'
    if(os.path.exists(dropbox_path + paper + 'summary_age_calibrators.tex')):
        os.remove(dropbox_path + paper + 'summary_age_calibrators.tex')
    file_sources = open(dropbox_path + paper + 'summary_age_calibrators.tex','x')

    #Header
    file_sources.write('\\begin{deluxetable*}{lccc}[ht!]\n')
    file_sources.write('\\tablewidth{290pt}\n')
    file_sources.write('\\tabletypesize{\scriptsize}\n')
    file_sources.write('\\tablecaption{Age Calibrators summary. \\label{table:age_cal}}\n')
    file_sources.write('\\tablehead{\\colhead{Reference} & \\colhead{Spectral} & \\colhead{Stars in Moving groups} & \\colhead{Stars with white dwarf companion} \\\ & Resolution & Compatible/All & Compatible/All\n}')
    file_sources.write('\\startdata \n')
    
    for x,y,z in zip(source_ref_table['source_ref'],source_ref_table['cite'],
                     source_ref_table['resolution']):
        mask_source = age_calibrators['source_ref'] == x
        mask_wd = age_calibrators['group_num'] == 0
        mask_mg = age_calibrators['group_num'] != 0
        n_wd = len(age_calibrators[mask_source*mask_wd*mask_ha])
        n_mg = len(age_calibrators[mask_source*mask_mg*mask_ha])
        n_wd_all = len(age_calibrators[mask_source*mask_wd])
        n_mg_all = len(age_calibrators[mask_source*mask_mg])
        if(any(np.array([n_wd,n_mg,n_wd_all,n_mg_all])!=0)):
            if(all(np.array([n_wd,n_wd_all])==0)):
                if(n_mg==n_mg_all):
                    file_sources.write(y+'&'+z+'&'+str(n_mg)+'&'+'-'+'\\\ \n')
                else:
                    file_sources.write(y+'&'+z+'&'+str(n_mg)+'/'+str(n_mg_all)+'&'
                                       +'-'+'\\\ \n')
            elif(n_mg==n_mg_all and n_wd==n_wd_all):
                file_sources.write(y+'&'+z+'&'+str(n_mg)+'&'
                                   +str(n_wd)+'\\\ \n') 
            else:
                file_sources.write(y+'&'+z+'&'+str(n_mg)+'/'+str(n_mg_all)+'&'
                                   +str(n_wd)+'/'+str(n_wd_all)+'\\\ \n') 

    file_sources.write('\\enddata \n')
    file_sources.write('\\end{deluxetable*}\n')

    return 0

#binaries = Table.read('../Catalogs/wdm_binaries.fits')
#make_table_for_wd(binaries)
#make_table_wd_ages(binaries)
    
#age_calibrators = Table.read('../Catalogs/age_calibrators_bayes.fits')
#make_table_summary_age_calibrators(age_calibrators)