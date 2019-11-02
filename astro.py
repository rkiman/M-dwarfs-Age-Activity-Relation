import numpy as np
from astropy import units as u

def calc_abs_mag(m, dist):
	'''
	This function calculates absolute magnitude from aparent magnitude.

	needs:
	m (array/float) aparent magnitude
	dist (array/float) distance (pc)
	'''
	mask = np.array([True if(str(x)!= 'nan' and str(x)!= 'NaN' and str(x)!= 'Inf' and str(x)!= 'inf' and x>0) else False for x in dist])

	abs_mag = np.zeros(len(dist))
	abs_mag[mask] = m[mask] + 5.0 - 5.0*np.log10(dist[mask])
	abs_mag[~mask] = np.nan

	return abs_mag
