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
    os.remove('data/wdm_co_movers_full.csv')
    wdm_table = Table()
    wdm_table['source_id'] = binaries['Source']
    wdm_table['TeffH'] = binaries['TeffH']
    wdm_table['e_TeffH'] = binaries['e_TeffH']
    wdm_table['loggH'] = binaries['loggH']
    wdm_table['e_loggH'] = binaries['e_loggH']
    wdm_table['ms_age_median_yr'] = binaries['ms_age_median']
    wdm_table['ms_age_err_low_yr'] = binaries['ms_age_err_low']
    wdm_table['ms_age_err_high_yr'] = binaries['ms_age_err_high']
    wdm_table['cooling_age_median_yr'] = binaries['cooling_age_median']
    wdm_table['cooling_age_err_low_yr'] = binaries['cooling_age_err_low']
    wdm_table['cooling_age_err_high_yr'] = binaries['cooling_age_err_high']
    wdm_table['total_age_median_yr'] = binaries['total_age_median']
    wdm_table['total_age_err_low_yr'] = binaries['total_age_err_low']
    wdm_table['total_age_err_high_yr'] = binaries['total_age_err_high']
    wdm_table['initial_mass_median'] = binaries['initial_mass_median']
    wdm_table['initial_mass_err_low'] = binaries['initial_mass_err_low']
    wdm_table['initial_mass_err_high'] = binaries['initial_mass_err_high']
    wdm_table['final_mass_median'] = binaries['final_mass_median']
    wdm_table['final_mass_err_low'] = binaries['final_mass_err_low']
    wdm_table['final_mass_err_high'] = binaries['final_mass_err_high']
    wdm_table.write('data/wdm_co_movers_full.csv',format='csv')

    os.remove('data/wdm_co_movers_ages.csv')
    wdm_table = Table()
    wdm_table['source_id_m'] = binaries['m_source_id']
    wdm_table['source_id_wd'] = binaries['Source']
    wdm_table['total_age_median_yr'] = binaries['total_age_median']
    wdm_table['total_age_err_low_yr'] = binaries['total_age_err_low']
    wdm_table['total_age_err_high_yr'] = binaries['total_age_err_high']
    wdm_table.write('data/wdm_co_movers_ages.csv',format='csv')
    
    dropbox_path = '/Users/rociokiman/Dropbox/Apps/Overleaf/'
    paper = 'wdwarfdate/'
    os.remove(dropbox_path + paper + 'wd_summary.tex')
    file_sources = open(dropbox_path + paper + 'wd_summary.tex','x')
    
    n_wd = 10
    rand_idx = np.random.randint(0,len(binaries),n_wd)

    #Header
    file_sources.write('\\begin{deluxetable*}{cccccccc}[ht!]\n')
    file_sources.write('\\tablewidth{290pt}\n')
    file_sources.write('\\tabletypesize{\scriptsize}\n')
    file_sources.write('\\tablecaption{WD sample. \\label{table:wd_sample}}\n')
    file_sources.write('\\tablehead{\
                        \\colhead{\\textit{Gaia} source id}\
                        &\\colhead{$T_{\\rm eff}$}\
                        & \\colhead{$\log (g)$}\
                        & \\colhead{$t_{\\rm ms}$(Gyr)}\
                        & \\colhead{$t_{\\rm cool}$(Gyr)}\
                        & \\colhead{$t_{\\rm tot}$(Gyr)}\
                        & \\colhead{$m_{\\rm i}$(M\\textsubscript{\(\odot\)})}\
                        & \\colhead{$m_{\\rm f}$(M\\textsubscript{\(\odot\)})}\
                        \n}')
    file_sources.write('\\startdata \n')
    
    for i in rand_idx:
        #source id
        source_id_i = str(binaries['Source'][i])
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
        file_sources.write(source_id_i +'&'+ teff+'&'+logg+'&'+ms+'&'+cool+'&'
                           +tot+'&'+mini+'&'+mfin+'\\\ \n') 

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
    
    n_wd = 10
    rand_idx = np.random.randint(0,len(binaries),n_wd)

    #Header
    file_sources.write('\\begin{deluxetable*}{ccc}[ht!]\n')
    file_sources.write('\\tablewidth{290pt}\n')
    file_sources.write('\\tabletypesize{\scriptsize}\n')
    file_sources.write('\\tablecaption{Ages for the white dwarfs co-moving with an M dwarf. \\label{table:wd_sample}}\n')
    file_sources.write('\\tablehead{\
                        \\multicolumn{2}{c}{\\textit{Gaia} source id} \
                        & \\colhead{Total age} \
                        \\\ M dwarf & White dwarf & (Gyr) \n}')
    file_sources.write('\\startdata \n')
    
    for i in rand_idx:
        source_id_m = '$' + str(np.round(binaries['m_source_id'][i],2)) + '$'
        source_id_wd = '$' + str(np.round(binaries['Source'][i],2)) + '$'
        age = str(np.round(binaries['total_age_median'][i]/1e9,2))
        err_low = str(np.round(binaries['total_age_err_low'][i]/1e9,2))
        err_high = str(np.round(binaries['total_age_err_high'][i]/1e9,2))
        tot = '$' + age + '_{-' + err_low + '}^{+' + err_high + '}$'
        file_sources.write(source_id_m+'&'+source_id_wd+'&'+tot+'\\\ \n') 

    file_sources.write('\\enddata \n')
    file_sources.write('\\end{deluxetable*}\n')
    
    data_path='/Users/rociokiman/Documents/M-dwarfs-Age-Activity-Relation/data'
    if(os.path.exists(data_path + '/wd_ages.csv')):
        os.remove(data_path + '/wd_ages.csv')
    data_table = open(data_path + '/wd_ages.csv','x')
    data_table.write('#Gaia id md\t Gaia id wd\t total age (yr)\t error low (yr)\t error high (yr)\n')
    for x,y,z,z1,z2 in zip(binaries['m_source_id'],binaries['Source'],
                           binaries['total_age_median'],
                           binaries['total_age_err_low'],
                           binaries['total_age_err_high']):
        data_table.write(str(x)+'\t'+str(y)+'\t'+str(z)+'\t'+str(z1)
                         +'\t'+str(z2)+'\n')
    data_table.close()
    
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
    
    n_tot_all = np.array([int(data_compatible[:,1][data_compatible[:,0] == x]) 
                          for x in source_ref_table['source_num']])
    mask_sort_n = np.flip(np.argsort(n_tot_all))
    
    for x,y,z in zip(source_ref_table['source_num'][mask_sort_n],
                     source_ref_table['cite'][mask_sort_n],
                     source_ref_table['resolution'][mask_sort_n]):
        if(x==25):
            y='LG11'+'\\tablenotemark{d}'
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
        if(x==19):
            n_comp=n_tot 
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
    file_sources.write('\\tablenotetext{d}{\citet{Lepine2013,Gaidos2014} with \
                       additional data observed in \
                       an identical manner.}\n')
    file_sources.write('\\end{deluxetable*}\n')

    return 0

def make_table_accretors(accretors):
    '''
    Makes a latex table with a summary of the accreators found in the
    literature search sample using the criterion from 
    White,R.J. & Basri,G. Astrophys. J. 582, 1109–1122 (2003).
    '''
    n_acc = 10
    
    if(os.path.exists('data/possible_accretors.csv')):
        os.remove('data/possible_accretors.csv')  
    table_acc = Table()
    table_acc['source_id'] = accretors['source_id']
    table_acc['spt'] = accretors['spt']
    table_acc['ewha'] = accretors['ewha_all']
    table_acc['ewha_error'] = accretors['ewha_error_all']
    d=[calc_delta_ha_for_accretors(x,y) for x,y in zip(table_acc['spt'],
                                                       table_acc['ewha'])]
    delta_ha = np.array([np.round(x,2) for x in d])
    table_acc['delta_ha'] = delta_ha
    table_acc['star_index'] = accretors['star_index']

    mask_delta = delta_ha >= 0
    assert all([x in table_acc['star_index'][mask_delta] for x in table_acc['star_index'][~mask_delta]])
    table_acc = table_acc[mask_delta]
    table_acc.write('data/possible_accretors.csv',format='csv')
        
    delta_ha = np.ones(n_acc)*-1
    spt = np.ones(n_acc)*-1
    ewha_error = np.ones(n_acc)*np.nan
    while((any(delta_ha<0) or any(np.isnan(ewha_error)) or any(spt<0))):
        rand_idx = np.random.randint(0,len(accretors),n_acc)    
        source_id = accretors['source_id'][rand_idx]
        spt = accretors['spt'][rand_idx]
        ewha = accretors['ewha_all'][rand_idx]
        ewha_error = accretors['ewha_error_all'][rand_idx]
        d = [calc_delta_ha_for_accretors(x,y) for x,y in zip(spt,ewha)]
        delta_ha = np.array([np.round(x,2) for x in d])
    
    
    dropbox_path = '/Users/rociokiman/Dropbox/Apps/Overleaf/'
    paper = 'Age-Activity Relation for M dwarfs/'
    if(os.path.exists(dropbox_path + paper + 'summary_accretors.tex')):
        os.remove(dropbox_path + paper + 'summary_accretors.tex')
    file_sources = open(dropbox_path + paper + 'summary_accretors.tex','x')
    
    #Header
    header_text = '\\tablehead{\
    \\colhead{\\textit{Gaia} source id} \
    & \\colhead{SpT} \
    & \\colhead{$\\haew$} \
    & \\colhead{$\\Delta \\haew$\\tablenotemark{a}} \n}'
    
    title = '\\tablecaption{Short sample of $\haew$ outliers,\
    possibly accreting according to the \
    criterion from \\citet{White2003}. The full sample can be found online.\
    \\label{table:accretors}}\n'
    
    file_sources.write('\\begin{deluxetable*}{lccc}[ht!]\n')
    file_sources.write('\\tablewidth{290pt}\n')
    file_sources.write('\\tabletypesize{\scriptsize}\n')
    file_sources.write(title)
    file_sources.write(header_text)
    file_sources.write('\\startdata \n')
    
    
    for x,s,y,z,w in zip(source_id,spt,ewha,ewha_error,delta_ha):
        file_sources.write(str(x) +' &'+
                           'M'+str(np.round(s,1))+' & $'+ 
                           str(np.round(y,2)) + '\\pm' 
                           +str(np.round(z,2)) +'$ & $'
                           +str(w) +'$ \\\ \n')
    file_sources.write('\\enddata \n')
    file_sources.write('\\tablenotetext{a}{Delta above the $\\haew$ limit.}\n')
    file_sources.write('\\end{deluxetable*}\n')
    
