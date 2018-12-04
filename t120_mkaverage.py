import glob
from astropy import units as u
from astropy.io import fits
from ccdproc import CCDData, Combiner
import numpy as np
from ccdproc import Combiner, ImageFileCollection,fits_ccddata_writer
from ccdproc import wcs_project
from astropy import wcs
import astroalign as aa
from t120_init import t120


def t120_mkaverage(imgdir=t120.t120_redu_dir,img_root_name='',img_ext_name='-c.fits',
                    master_root_name='master'):

    dict_ccd_data	= {}
    dict_ccd_proj	= {}
    
    # get list of images
    list_img		= glob.glob(imgdir+img_root_name+'*'+img_ext_name)
    
    # load first image for reference
    hdulist		= fits.open(list_img[0])
    ref_img_array	= np.array(hdulist[0].data)
    
    # now loops over all the files
    for fit_file in list_img:
     t120.log.info('now opening file: '+fit_file)
     hdulist = fits.open(fit_file)
     filter_name = hdulist[0].header['FILTER']
     t120.log.info('filter='+filter_name)
     if not dict_ccd_data.has_key(filter_name):
        dict_ccd_proj[filter_name] = []
        dict_ccd_proj[filter_name].append(CCDData(hdulist[0].data,unit=u.adu))
     else:
         proj_img_array = aa.register(np.array(hdulist[0].data),ref_img_array)
         dict_ccd_proj[filter_name].append(CCDData(proj_img_array,unit=u.adu))
    
    t120.log.info('now loop over the filters')
    for filter_name in dict_ccd_proj:
     combiner_proj = Combiner(dict_ccd_proj[filter_name])     
     master_proj = combiner_proj.median_combine()
     hdu_proj = master_proj.to_hdu()
     master_file_proj = imgdir+img_root_name+'-'+master_root_name+'-'+filter_name+'-proj.fits'
     hdu_proj.writeto(master_file_proj,overwrite=True)
     t120.log.info('Master img saved in '+master_file_proj)
    return

#t120_mkaverage(imgdir='/Volumes/BackUpMac/OHP/2017/2017-08-21/Astrometry/',img_root_name='m33')
t120_mkaverage(imgdir='/Volumes/BackUpMac/OHP/2017/2017-08-22/Astrometry/',img_root_name='M74')
