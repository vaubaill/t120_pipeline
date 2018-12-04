"""Clean Astrometry directory from file for which the SCAMP contrast is bad.

"""
import os
from astropy.io.votable import parse
from t120_init import t120

# user defined
mincontrast=5.0


scampxml = t120.t120_scampres_dir + 'scamp.xml'

if not os.path.isfile(scampxml):
    msg = '*** FATAL ERROR: file '+scampxml+' does not exist.'
    t120.log.error(msg)
    raise IOError(msg)
    
votable = parse(scampxml)
table = votable.get_first_table()
contrast = (table.array['XY_Contrast']).data
catname = (table.array['Catalog_Name']).data

for (contrast,catname) in zip((table.array['XY_Contrast']).data,(table.array['Catalog_Name']).data):
    if (contrast < mincontrast):
        print('File '+catname+' has a bad contrast of '+str(contrast)+' please delete file: '+t120.t120_astr_dir +catname.replace('ldac','fits'))
        fitsfile = t120.t120_astr_dir + catname.replace('ldac','fits')
        if os.path.exists(fitsfile):
            t120.log.info('Now deleting file: '+fitsfile)
            os.remove(fitsfile)
        else:
            msg = 'Impossible to delete file: '+fitsfile
            t120.log.warning(msg)
            raise IOError(msg)

