import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

def WISE_id_to_ra_dec(wise_name):
    ra_all,dec_all = [],[]
    for x in wise_name:
        
        x.replace('J', '')
        
        ra_name='none'
        dec_name='none'
        
        if('A' in x or 'B' in x):
            ra_all.append(np.nan)
            dec_all.append(np.nan)
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
            
            ra_all.append(coord.ra.deg)
            dec_all.append(coord.dec.deg)
            
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
            
            ra_all.append(coord.ra.deg)
            dec_all.append(coord.dec.deg)
            
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