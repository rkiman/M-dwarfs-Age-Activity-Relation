import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from scipy.interpolate import interp1d

def WISE_id_to_ra_dec(wise_name):
    N = len(wise_name)
    ra_all,dec_all = np.ones(N)*np.nan, np.ones(N)*np.nan
    
    for x,i in zip(wise_name,range(N)):

        x = x.replace('J', '')
        x = x.replace('.', '')

        ra_name='none'
        dec_name='none'
        
        if('A' in x or 'B' in x):
            ra_all[i] = np.nan
            dec_all[i] = np.nan
        elif('+' in x):
            ra,dec = x.split('+')

            if(len(ra)==8 or len(ra)==7):
                ra_name = ra[:2]+':'+ra[2:4]+':'+ra[4:6]+'.'+ra[6:]
            elif(len(ra)==6):
                ra_name = ra[:2]+':'+ra[2:4]+':'+ra[4:6]
            else:
                print(ra)
                
            if(len(dec)==8 or len(dec)==7):
                dec_name = dec[:2]+':'+dec[2:4]+':'+dec[4:6]+'.'+dec[6:]
            elif(len(dec)==6):
                dec_name = dec[:2]+':'+dec[2:4]+':'+dec[4:6]
            else:
                print(dec)
            star_name = ra_name+'+'+dec_name
            #print(star_name)
            coord = SkyCoord(star_name, unit=(u.hourangle, u.deg))
            
            ra_all[i] = coord.ra.deg
            dec_all[i] = coord.dec.deg
            
        elif('-' in x):
            ra,dec = x.split('-')
            if(len(ra)==8 or len(ra)==7):
                ra_name = ra[:2]+':'+ra[2:4]+':'+ra[4:6]+'.'+ra[6:]
            elif(len(ra)==6):
                ra_name = ra[:2]+':'+ra[2:4]+':'+ra[4:6]
            else:
                print(ra)
                
            if(len(dec)==8 or len(dec)==7):
                dec_name = dec[:2]+':'+dec[2:4]+':'+dec[4:6]+'.'+dec[6:]
            elif(len(dec)==6):
                dec_name = dec[:2]+':'+dec[2:4]+':'+dec[4:6]
            else:
                print(dec)
            star_name = ra_name+'-'+dec_name
            #print(star_name)
            coord = SkyCoord(star_name, unit=(u.hourangle, u.deg))
            
            ra_all[i] = coord.ra.deg
            dec_all[i] = coord.dec.deg
            
    return np.array(ra_all), np.array(dec_all)


def calc_abs_mag(m, dist):
	'''
	This function calculates absolute magnitude from aparent magnitude.

	needs:
	m (array/float) aparent magnitude
	dist (array/float) distance (pc)
	'''
	mask = np.array([True if(str(x)!='nan' 
                          and str(x)!='NaN' 
                          and str(x)!='Inf' 
                          and str(x)!='inf' and x>0) else False for x in dist])

	abs_mag = np.zeros(len(dist))
	abs_mag[mask] = m[mask] + 5.0 - 5.0*np.log10(dist[mask])
	abs_mag[~mask] = np.nan

	return abs_mag



def hstodeg(ra_hs,dec_hs):
    ra,dec = [],[]
    for i in range(len(ra_hs)):
        x,y = ra_hs[i],dec_hs[i]
        if(str(x)=='--' or str(y)=='--'):
            ra.append(np.nan)
            dec.append(np.nan)        
        else:
            radec_hs_i = x + ' ' +  y
            c = SkyCoord(radec_hs_i, unit=(u.hourangle, u.deg))
            ra.append(c.ra.deg)
            dec.append(c.dec.deg)
    return np.array(ra),np.array(dec)

def calc_lhalbol(ewha,ewha_error,spt):
    chi_douglas2014 = np.array([6.6453, 6.0334, 5.2658, 4.4872, 3.5926, 
                                2.4768, 1.7363, 1.2057, 0.6122, 0.3522])*1e-5
    spt_douglas2014 = np.arange(0,10)
    f_chi = interp1d(spt_douglas2014,chi_douglas2014)

    N = len(ewha)
    lhalbol = np.ones(N)*np.nan 
    lhalbol_err = np.ones(N)*np.nan
    for i in range(N):
        if((0 <= spt[i]) and (spt[i] <= 9)):
            if(~np.isnan(ewha_error[i]+ewha[i])):
                dist_ewha = np.random.normal(ewha[i],ewha_error[i],2000)
                dist_lhalbol = dist_ewha*f_chi(spt[i])
                lhalbol[i] = np.nanmedian(dist_lhalbol)
                lhalbol_err[i] = np.nanstd(dist_lhalbol)
            elif(np.isnan(ewha_error[i]) and ~np.isnan(ewha[i])):
                lhalbol[i] = ewha[i]*f_chi(spt[i])
    
    return lhalbol,lhalbol_err

def calc_log_lhalbol(ewha,ewha_error,spt):
    chi_douglas2014 = np.array([6.6453, 6.0334, 5.2658, 4.4872, 3.5926, 
                                2.4768, 1.7363, 1.2057, 0.6122, 0.3522])*1e-5
    spt_douglas2014 = np.arange(0,10)
    f_chi = interp1d(spt_douglas2014,chi_douglas2014)

    N = len(ewha)
    log_lhalbol = np.ones(N)*np.nan 
    log_lhalbol_err = np.ones(N)*np.nan
    for i in range(N):
        if((0 <= spt[i]) and (spt[i] <= 9)):
            if(~np.isnan(ewha_error[i]+ewha[i])):
                dist_ewha = np.random.normal(ewha[i],ewha_error[i],2000)
                dist_log_lhalbol = np.log10(dist_ewha*f_chi(spt[i]))
                log_lhalbol[i] = np.nanmedian(dist_log_lhalbol)
                log_lhalbol_err[i] = np.nanstd(dist_log_lhalbol)
            elif(np.isnan(ewha_error[i]) and ~np.isnan(ewha[i])):
                log_lhalbol[i] = np.log10(ewha[i]*f_chi(spt[i]))
    
    return log_lhalbol,log_lhalbol_err
            
def organize_table_format(columns):
    
    labels = ['ra', 'dec', 'spt', 'gaia_source_id', 'ra_gaia', 'dec_gaia', 
              'pmra', 'pmra_error', 'pmdec', 'pmdec_error', 'parallax', 
              'parallax_error', 'phot_g_mean_flux', 'phot_g_mean_flux_error',
              'phot_g_mean_mag', 'phot_rp_mean_flux', 
              'phot_rp_mean_flux_error', 'phot_rp_mean_mag', 
              'phot_bp_mean_flux', 'phot_bp_mean_flux_error', 
              'phot_bp_mean_mag', 'g_corr', 'rp_corr', 'ewha', 
              'ewha_error', 'ewha_all', 
              'ewha_error_all', 'lhalbol', 'lhalbol_error', 'age', 
              'age_error_low', 'age_error_high', 'group_num', 'group_name',
              'star_index',
              'source_num', 'source_ref']
    N = len(labels)
    
    organized_table = Table()
    for i in range(N):
        organized_table[labels[i]] = np.array(columns[i])

    return organized_table


def add_corrected_magnitudes(catalog):
    deredden = Table.read('Catalogs/Sources/deredden_gaia_mags.fits')
    N = len(catalog['source_id'])
    g_corr = np.ones(N)*np.nan
    rp_corr = np.ones(N)*np.nan
    
    for i in range(N):
        mask = catalog['source_id'][i] == deredden['SOURCE_ID']
        if(len(deredden['G_CORR'][mask])!=0):
            g_corr[i] = deredden['G_CORR'][mask][0]
            rp_corr[i] = deredden['GRP_CORR'][mask][0]
    
    catalog['g_corr'] = g_corr
    catalog['rp_corr'] = rp_corr
    
    return catalog

def calc_number_single_stars(catalog):
    '''
    We decided to use all the measurements we had for ha. But we need to 
    report how many single stars we have in each case. This is what this
    function calculates using the star_index column. If two stars have the 
    same number in that column then they are the same star.
    '''
    star_index = catalog['star_index']
    n_single = len(set(star_index))
    
    return n_single

def find_idx_singles(catalog):
    same_star = catalog['star_index']
    seen = []
    idx = []
    for i,x in enumerate(same_star):
        if(x=='nan'):
            idx.append(i)
        else:
            if(x not in seen):
                seen.append(x)
                idx.append(i)
    return np.array(idx)

def color_to_spt(color):
    spt_model = np.linspace(-2,14,100) 
    g_rp_model = -0.0036*spt_model**2+0.11*spt_model+0.89
    f = interp1d(g_rp_model, spt_model)
    return f(color)

def get_spt(spt,g_rp):
    mask_color = (g_rp > 0.6556) * (g_rp < 1.7244000000000002) * np.isnan(spt)
    spt_from_color_1 = color_to_spt(g_rp[mask_color])
    spt_from_color = np.ones(len(g_rp))*np.nan
    spt_from_color[mask_color] = spt_from_color_1
    spt_from_color[~np.isnan(spt)] =spt[~np.isnan(spt)]
    return spt_from_color
