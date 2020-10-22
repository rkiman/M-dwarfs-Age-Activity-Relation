#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.table import Table
import numpy as np
import os
from astropy.coordinates import SkyCoord
import astropy.units as u
from datetime import datetime
from .astro import calc_number_single_stars

def find_new_members():
    
    #Open log file to record result numbers
    log_file = open('log.txt','a')
    log_file.write('\n')
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    log_file.write("Log's date and time: {}\n".format(dt_string))
    log_file.flush()
    
    #Load table with results from banyan
    all_groups = Table.read('Catalogs/literature_search_all_groups.fits')
    mask_not_mem = (all_groups['group_name']!='nan') * (all_groups['good_mem']==0)
    n_rejected=calc_number_single_stars(all_groups[mask_not_mem])
    text='Number of single stars that were rejected as members with banyan: {}\n'
    log_file.write(text.format(n_rejected))
    
    #Good members were identified with the column good_mem == 1
    #Also making the sums with the compatible members
    mask_mem = (all_groups['good_mem'] == 1) * (~np.isnan(all_groups['ewha']))
    
    #make a sub table with good members only
    all_mem = all_groups[mask_mem]
    #Columns I'm going to need to select new members
    original_group = all_mem['group_name']
    banyan_group = all_mem['best_ya']
    star_index = all_mem['star_index']
    
    new_mem = []
    seen = []
    change_mem = 0
    same_mem = 0
    #Look over the number star_index which indicates which stars are repeated
    for i in range(len(all_mem)):
        x = star_index[i]
        #star_index = nan indicates that the star is not repeated
        if(np.isnan(x)):
            #Compare original classification to banyan's classification
            if(original_group[i]!=banyan_group[i]):
                #original_group = nan indicates that the star originally 
                #did not have a group
                if(original_group[i]=='nan'):
                    new_mem.append(i)
                else:
                    change_mem+=1
                    new_mem.append(i)
            elif(original_group[i]==banyan_group[i]):
                same_mem +=1
        else:
            #If the star is repeated I check all of the repeated measurements 
            #at the same time
            if(x not in seen):
                seen.append(x)
                mask_star_index = star_index == x
                if(len(star_index[mask_star_index])>0): 
                    original_group_ss = original_group[mask_star_index]
                    banyan_group_ss = banyan_group[mask_star_index]
                    if(any(original_group_ss!=banyan_group_ss)):
                        #If none of the repeated stars have a classification 
                        #then is new
                        if(all(original_group_ss=='nan')):
                            #Check that all banyan's classifications agree for 
                            #the same star
                            assert all(banyan_group_ss==banyan_group_ss[0])
                            new_mem.append(i)
                        #If any of the repeated stars had a classification, 
                        #then it means the
                        #Banyan classification changed the group.
                        else:
                            assert all(banyan_group_ss==banyan_group_ss[0])
                            groups_i = original_group_ss[original_group_ss!='nan']
                            if(any(groups_i==banyan_group_ss[0])):
                                same_mem +=1
                            else:
                                change_mem+=1
                                new_mem.append(i)
                            #The rest kept the classification they had
                            #else:
                            #    print(x)
                            #    for w,y in zip(original_group_ss,banyan_group_ss):
                            #        print(w,y)
                    elif(all(original_group_ss==banyan_group_ss)):
                        same_mem +=1
                            
    #Check that the members counted are single
    assert len(set(star_index[new_mem])) == len(new_mem)
    
    text = 'Number of single members that changed their group with banyan: {}\n'
    log_file.write(text.format(change_mem))
    text = 'Number of single members that did not change their group: {}\n'
    log_file.write(text.format(same_mem))
    text = 'Number of single members that are "new" members: {}\n'
    log_file.write(text.format(len(new_mem)-change_mem))
    
    
    #Check if the new members are in the literature
    #Load columns I need only for the new members according to my literature search
    #sample
    ra_mg = all_mem['ra'][new_mem]
    dec_mg = all_mem['dec'][new_mem]
    group_mg = all_mem['best_ya'][new_mem]
    source_mg = all_mem['source_ref'][new_mem]
    source_id_mg = all_mem['source_id'][new_mem]
    spt_mg = all_mem['spt'][new_mem]
    
    #Load literature catalogs with young associations members to check if the
    #ones I found are new or not
    BF_mem = Table.read('Catalogs/Sources/BF_members.fit')
    all_mem_banyan = Table.read('Catalogs/Sources/allmembers_banyan.fits')
    new_gf = Table.read('Catalogs/Sources/GF2018.fit')
    roeser2011 = Table.read('Catalogs/Sources/Roeser2011.fit')
    
    ra_roeser2011 = np.array(roeser2011['RAJ2000'])
    dec_roeser2011 = np.array(roeser2011['DEJ2000'])
    c_roeser2011 = SkyCoord(ra=ra_roeser2011*u.deg, dec=dec_roeser2011*u.deg)
    
    ra_bf_mem = np.array(BF_mem['RAJ2000'])
    dec_bf_mem = np.array(BF_mem['DEJ2000'])
    c_bf_mem = SkyCoord(ra=ra_bf_mem*u.deg, dec=dec_bf_mem*u.deg)
    
    ra_all_mem = all_mem_banyan['ra']
    dec_all_mem = all_mem_banyan['dec']
    c_all_mem = SkyCoord(ra=ra_all_mem*u.deg, dec=dec_all_mem*u.deg)
    
    ra_new_gf = np.array(new_gf['RAJ2000'])
    dec_new_gf = np.array(new_gf['DEJ2000'])
    c_new_gf = SkyCoord(ra=ra_new_gf*u.deg, dec=dec_new_gf*u.deg)
    
    #File to save data from new members
    if(os.path.exists('data/new_members_data.csv')):
        os.remove('data/new_members_data.csv')
    new_members = open('data/new_members_data.csv','a')
    new_members.write('source_id\tra\tdec\tspt\tgroup\treference\n')
    
    tot_match = 0
    tot_new_mem = 0
    for ra_i,dec_i,group_i,source_i,source_id_i,spt_i in zip(ra_mg,dec_mg,group_mg,
                                                             source_mg,source_id_mg,
                                                             spt_mg):
        c1 = SkyCoord(ra=ra_i*u.deg, dec=dec_i*u.deg)
        separation_bf_mem = c1.separation(c_bf_mem).arcsec
        separation_all_mem = c1.separation(c_all_mem).arcsec
        separation_new_gf = c1.separation(c_new_gf).arcsec
        separation_roeser2011 = c1.separation(c_roeser2011).arcsec
        if(any(np.concatenate((separation_bf_mem,separation_all_mem,
                               separation_new_gf,separation_roeser2011))<=6)):
            tot_match+=1
        else:
            tot_new_mem+=1
            new_members.write(str(source_id_i)+'\t'+str(c1.ra.deg)
                              +'\t'+str(c1.dec.deg)+'\t'+str(spt_i)
                              +'\t'+str(group_i) 
                              + '\t'+str(source_i) + '\n')
            
    new_members.close()
    
    text = 'Number of single members that were identified by gagne and roeser: {}\n'
    log_file.write(text.format(tot_match))
    
    text = 'Number of single members that are posible new members: {}\n'
    log_file.write(text.format(tot_new_mem))
    
    log_file.close()
    
    
def make_table_new_members(dropbox_path,paper):
    
    #Make table for paper
    source_ref_table = Table.read('data/source_ref.csv')
    new_mem_table = Table.read('data/new_members_data.csv',format='csv',
                               delimiter='\t')

    if(os.path.exists(dropbox_path + paper + 'summary_new_mems.tex')):
        os.remove(dropbox_path + paper + 'summary_new_mems.tex')
    file_sources = open(dropbox_path + paper + 'summary_new_mems.tex','x')
    
    #Header
    header_text = '\\tablehead{\
    \\colhead{\\textit{Gaia} source id} \
    & \\colhead{SpT}\
    & \\colhead{group\\tablenotemark{a}}\
    & \\colhead{Reference} \n}'
    
    title = '\\tablecaption{Sample of new members. We show only $10$ stars as an example \\label{table:newmem}}\n'
    
    file_sources.write('\\begin{deluxetable*}{lccc}[ht!]\n')
    file_sources.write('\\tablewidth{290pt}\n')
    file_sources.write('\\tabletypesize{\scriptsize}\n')
    file_sources.write(title)
    file_sources.write(header_text)
    file_sources.write('\\startdata \n')
    
    n_mem = 10
    rand_idx = np.random.randint(0,len(new_mem_table),n_mem)
    source_id = new_mem_table['source_id'][rand_idx]
    spt = new_mem_table['spt'][rand_idx]
    group = new_mem_table['group'][rand_idx]
    ref = new_mem_table['reference'][rand_idx]
    reference_source = source_ref_table['source_ref']
    cite = [source_ref_table['cite'][reference_source==x][0] for x in ref]
    
    for x,y,z,w in zip(source_id,spt,group,cite):
        file_sources.write(str(x)+'& M$'+str(np.round(y,1))+'$&'+str(z)
        +'&'+str(w)+' \\\ \n')
    file_sources.write('\\enddata \n')
    #note = '\\tablenotetext{a}{Total number of new members per group:'
    #for x,y in zip(groups_new,n_groups):
    #    note = note + '\t' + str(x) + ':' + str(y) +', '
    #file_sources.write(note+ '.}\n')
    file_sources.write('\\end{deluxetable*}\n')
    file_sources.close()
    
def make_table_summary_members(dropbox_path,paper):
    #Summary of members
    mg_confirmed = Table.read('Catalogs/literature_search_mg.fits')
    group_num = mg_confirmed['group_num']
    ewha = mg_confirmed['ewha']
    new_mem_table = Table.read('data/new_members_data.csv',format='csv',
                               delimiter='\t')
    groups_new = set(new_mem_table['group'])
    n_groups = [len(new_mem_table[new_mem_table['group']==x]) for x in groups_new]
    
    mg_ref = Table.read('data/moving_groups_ref.csv')
    mg_ref = mg_ref[:-1]
    idx_age = np.argsort(mg_ref['age'])
    
    if(os.path.exists(dropbox_path + paper + 'summary_all_mems.tex')):
        os.remove(dropbox_path + paper + 'summary_all_mems.tex')
    file_sources = open(dropbox_path + paper + 'summary_all_mems.tex','x')
    
    #Header
    header_text = '\\tablehead{\
    \\colhead{Young Association} \
    & \\colhead{Total members}\
    & \\colhead{New members} \n}'
    
    title = '\\tablecaption{Summary of M dwarfs in young associations. \\label{table:summary_mem}}\n'
    
    file_sources.write('\\begin{deluxetable*}{lcc}[ht!]\n')
    file_sources.write('\\tablewidth{290pt}\n')
    file_sources.write('\\tabletypesize{\scriptsize}\n')
    file_sources.write(title)
    file_sources.write(header_text)
    file_sources.write('\\startdata \n')
    
    n_groups = np.array(n_groups)
    for name_group_i,num_group_i in zip(mg_ref['name'][idx_age],
                                        mg_ref['group_num'][idx_age]):
        mask = (group_num == num_group_i) * (~np.isnan(ewha))
        n_group_tot = calc_number_single_stars(mg_confirmed[mask])
        n_group_new = 0
        #print(name_group_i)
        if(name_group_i in groups_new):
            mask_for_new_mem = np.array([x==name_group_i for x in groups_new])
            n_group_new = n_groups[mask_for_new_mem][0]
        if(n_group_tot!=0):
            file_sources.write(name_group_i+'& $'+str(n_group_tot)+'$&$'+str(n_group_new) +'$ \\\ \n')
            
    file_sources.write('\\enddata \n')
    file_sources.write('\\end{deluxetable*}\n')
    file_sources.close()