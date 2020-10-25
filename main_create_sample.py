
from astropy.table import Table
import os
import src 

#Identify compatible sample
file = 'Catalogs/literature_search_gaia_compatible.fits'
if(os.path.exists(file)):
    literature_search_gaia_compatible = Table.read(file)
else:
    literature_search_gaia_compatible = src.compile_compatible_sample()

#Compile sample of calibrators
file = 'Catalogs/age_calibrators.fits'
if(os.path.exists(file)):
    age_calibrators = Table.read(file)
else:
    age_calibrators = src.compile_age_calibrators(literature_search_gaia_compatible)

#identify photometric binaries
age_calibrators = src.add_binary_column(age_calibrators)

#Save sample
print('Saving age calibrators sample')
age_calibrators.write('Catalogs/age_calibrators.fits', overwrite=True)
print('Done sample of age calibrators')

#Identify new members of moving groups
file = 'data/new_members_data.csv'
if(~os.path.exists(file)):
    src.find_new_members()

#Make tables for the paper
dropbox_path = '/Users/rociokiman/Dropbox/Apps/Overleaf/'
paper = 'Age-Activity Relation for M dwarfs/'

binaries = Table.read('Catalogs/wdm_binaries.fits')
src.make_table_wd_ages(binaries)
    
age_calibrators = Table.read('Catalogs/age_calibrators.fits')
src.make_table_summary_age_calibrators(age_calibrators)
    
accretors = Table.read('Catalogs/literature_search_accretors.fits')
src.make_table_accretors(accretors)

src.make_table_new_members(dropbox_path,paper)

src.make_table_summary_members(dropbox_path,paper)