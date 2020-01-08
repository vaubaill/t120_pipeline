import os
import glob
from astropy.io.votable import parse

from t120_init import t120

def cleanastrodir(astrodir,mincontrast=5.0):
    """Remove all files for which the SCAMP contrast is not high enough.
    
    Parameters
    ----------
    astrodir : string
        root of the directory name where data are located.
    mincontrast : float, optional
        minimum contrast below which fits files are removed.
    
    Returns
    -------
    None.
    
    """
    scampxml = astrodir+t120.t120_scampres_dir + 'scamp.xml'
    
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
            t120.log.warning('File '+catname+' has a bad contrast of '+str(contrast)+' please delete file: '+t120.t120_astr_dir +catname.replace('ldac','fits'))
            fitsfile = t120.t120_astr_dir + catname.replace('ldac','fits')
            if os.path.exists(fitsfile):
                t120.log.info('Now deleting file: '+fitsfile)
                os.remove(fitsfile)
            else:
                msg = 'Impossible to delete file: '+fitsfile
                t120.log.warning(msg)
                raise IOError(msg)
    return

def cleancosmeticdir(rootdir):
    """Clean cosmetic directories from master files.
    
    Parameters
    ----------
    rootdir : string
        root directory name
    
    Returns
    -------
    None.
    
    """
    for subdir in [t120.t120_ofst_dir,t120.t120_flat_dir,t120.t120_dark_dir]:
        tocleandir = rootdir + t120.t120_ofst_dir
        t120.log.info('tocleandir: '+tocleandir)
        pattern = tocleandir+'master*.fit*'
        tocleanlist = glob.glob(pattern)
        if len(tocleanlist):
            for tocleanfile in tocleanlist:
                t120.log.info('Now removing file: '+tocleanfile)
                os.remove(tocleanfile)
    return
