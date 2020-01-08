import os
import glob
from astropy import units as u
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from ccdproc import CCDData, Combiner, subtract_dark,flat_correct,ImageFileCollection,ccd_process,fits_ccddata_reader,fits_ccddata_writer
import numpy as np
import ccdproc

from t120_init import t120

def t120_insert_radec(work_dir,ahead_file=t120.t120_scamp_ahead,reduc_root='-c.fits'):
    """Put object RA,DEC in fits header
    
    Parameters
    ----------
    work_dir : string
        directory name root where original fits image files are to be found.
    ahead_file : string
        Scamp ahead file to consider to put basic astrometry info.
    reduc_root : string
        Reduced file name suffix. Default is '-c.fit'.
    
    Returns
    -------
    None.
    
    
    """
    redu_dir=work_dir + t120.t120_redu_dir
    t120.log.info('redu_dir '+redu_dir)
    
    comm = {'OBJCTRA':'Nominal Right Ascension of center of image',
           'OBJCTDEC':'Nominal Declination of center of image'}
    
    for fit_file in glob.glob(redu_dir+'*-c.fits'):
        t120.log.info('*** Now opening image '+fit_file)
        hdulist = fits.open(fit_file)
        hdu = hdulist[0]
        hdr = hdu.header
        data = hdu.data
        target_name = hdr['OBJECT'].strip()
        hdr['OBJECT'] = (target_name.upper().replace(' ',''),'Target name')
        skycoo = SkyCoord(hdr['OBJCTRA'],hdr['OBJCTDEC'], unit=(u.hourangle, u.deg))
        time = Time(hdr['DATE-OBS'])
        site = EarthLocation(lat=hdr['SITELAT'], lon=hdr['SITELONG']) 
        hdr['CRVAL1'] = (skycoo.ra.to('deg').value ,'Reference Right ascencion in decimal deg')
        hdr['CRVAL2'] = (skycoo.dec.to('deg').value,'Reference Declination in decimal deg')
        # now put astrometry info from ahead_file
        aheader = fits.Header.fromfile(ahead_file,sep='\n',endcard=False,padding=False)
        for ahdr in aheader:
            #nhdr=len(hdr)
            t120.log.debug('Setting: '+ahdr+' with value: '+str(aheader.cards[ahdr].value))
            hdr.set(ahdr,value=aheader.cards[ahdr].value,comment=aheader.cards[ahdr].comment)
        hdr.add_history('BASIC ASTROMETRY WAS TAKEN FROM FILE '+ahead_file)
        hdu.writeto(fit_file,overwrite=True)
        t120.log.info('Image saved in '+fit_file)
    return

