import os
import glob
from astropy import units as u
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from ccdproc import CCDData, Combiner, subtract_dark,flat_correct,ImageFileCollection,ccd_process,fits_ccddata_reader,fits_ccddata_writer
import numpy as np
import ccdproc


def t120_insert_radec(work_dir,objcoo):

   comm		= {'OBJCTRA':'Nominal Right Ascension of center of image',
		   'OBJCTDEC':'Nominal Declination of center of image'}

   for fit_file in glob.glob(work_dir+'*-c.fits'):
    hdulist	= fits.open(fit_file)
    hdu		= hdulist[0]
    data	= hdu.data
    target_name	= hdu.header['OBJECT'].strip()
    hdu.header['OBJECT']		= (target_name.upper().replace(' ',''),'Target name')
    skycoo		= SkyCoord(objcoo[hdu.header['OBJECT']][0]+' '+objcoo[hdu.header['OBJECT']][1], unit=(u.hourangle, u.deg))
    time		= Time(hdu.header['DATE-OBS'])
    site		= EarthLocation(lat=hdu.header['SITELAT'], lon=hdu.header['SITELONG']) 
    altazframe		= AltAz(obstime=time, location=site)
    altaz		= skycoo.transform_to(altazframe)
    try:
      hdu.header['OBJCTRA']	= (objcoo[hdu.header['OBJECT']][0],comm['OBJCTRA'])
      hdu.header['OBJCTDEC']	= (objcoo[hdu.header['OBJECT']][1],comm['OBJCTDEC'])
      hdu.header['CRVAL1']		= (skycoo.ra.to('deg').value ,'Reference Right ascencion in decimal deg')
      hdu.header['CRVAL2']		= (skycoo.dec.to('deg').value,'Reference Declination in decimal deg')
      hdu.header['OBJCTALT']	= (altaz.alt.to('deg').value ,'Reference altitude (elevation) in deg')
      hdu.header['OBJCTAZ']	= (altaz.az.to('deg').value,'Reference azimuth in deg')
      hdu.header['OBJCTHA']	= ((altaz.az-180.0*u.deg).to('hourangle').value,'Hour angle')
      hdu.header['AIRMASS']	= (altaz.secz.value,'Air mass')
    except:
      msg= '*** FATAL ERROR: target '+target_name+' not found'
      raise ValueError(msg)
    hdu.writeto(fit_file,overwrite=True)


