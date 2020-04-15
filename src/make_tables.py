#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from astropy.table import Table
#import numpy as np

def make_table_for_paper_sources(source_num,Ncomp,compatible,total_comp,order):
     
    dropbox_path = '/Users/rociokiman/Dropbox/Apps/Overleaf/'
    paper = 'Age-Activity Relation for M dwarfs/'
    os.remove(dropbox_path + paper + 'sources_summary.tex')
    file_sources = open(dropbox_path + paper + 'sources_summary.tex','x')
    
    source_ref_table = Table.read('data/source_ref.csv')
    
    #Ntot = len(ra)
    #Ncomp = len(ls_compatible[~np.isnan(ls_compatible['ewha'])])
    n_sources = len(source_ref_table['source_ref'])
    total_incompatible = 0
    papers_notes = 'Other works checked in literature search but were considered incompatible: '
    #Header
    file_sources.write('\\begin{deluxetable*}{ccccc}[ht!]\n')
    file_sources.write('\\tablewidth{290pt}\n')
    file_sources.write('\\tabletypesize{\scriptsize}\n')
    file_sources.write('\\tablecaption{Compatible references for $\haew$. \\label{table:source_ref}}\n')
    file_sources.write('\\tablehead{\\colhead{Reference \\tablenotemark{a}} & \\colhead{Resolution} & \\colhead{$N$ Stars} & \\colhead{$N$ Stars Compatible} & \\colhead{OC \\tablenotemark{b}} \\\ & & Total='+str(Ncomp) +'& &\n}')
    file_sources.write('\\startdata \n')
    
    for i in range(n_sources):
        if i in compatible:
            mask = i == compatible
            compatible_i = int(total_comp[mask][0])
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
    file_sources.write('\\tablenotetext{b}{Order of compatibility. Order 1 is compatible with \\citet{Kiman2019}. Order 2 are compatible with order 1 catalogs.}\n')
    file_sources.write('\\end{deluxetable*}\n')

    return 0

