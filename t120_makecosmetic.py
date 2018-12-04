import os, sys, glob
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from ccdproc import CCDData, Combiner, subtract_overscan,ImageFileCollection,ccd_process,fits_ccddata_reader,fits_ccddata_writer
import numpy as np
import ccdproc
import silentremove
from t120_init import t120

def read_scamp_ahead(scampahead_file):
    try:
        with open(scampahead_file, 'r') as f:
            lines = f.read()
    except IOError as err:
        msg = '*** FATAL ERROR: file '+scampahead_file+' might not exist'
        t120.log.error(msg)
        raise IOError(msg)
    else:
        header = []
        for line in lines.split('\n'):
            if len(line):
                header.append(line.strip())
        del header[-1]	# remove 'END'
        return header


def put_scamp_header(header,scamp_header):
    for h in scamp_header:
        keyword = h.split('=')[0]
    try:
        value = float(h.split('=')[1].split('/')[0])
    except:
        value = h.split('=')[1].split('/')[0].replace('"', '').replace('\'', '')
    comm = h.split('=')[1].split('/')[1].strip('\'')
    header[keyword] = (value,comm)
    return header


def t120_mkoffset(offset_dir=t120.t120_ofst_dir,master_file_name=t120.t120_master_name):
    master_file = offset_dir+master_file_name
    listimg = ImageFileCollection(offset_dir) #,glob_include='*.fit',glob_exclude='*.fits')
    listccd = []
    for ccd,file_name in listimg.ccds(ccd_kwargs={'unit':'adu'},return_fname=True):
        t120.log.info('now considering file '+file_name)
        listccd.append(ccd)
    
    combiner = Combiner(listccd)
    t120.log.info('now making the combination')
    master_offset = combiner.median_combine()
    fits_ccddata_writer(master_offset,master_file)
    t120.log.info('Result saved in '+master_file)
    return

def t120_mkdark(flat_dir=t120.t120_flat_dir,master_file_name=t120.t120_master_name):
    master_file = flat_dir+master_file_name
    listimg = ImageFileCollection(flat_dir) #,glob_include='*.fit',glob_exclude='*.fits')
    listccd = []
    for ccd,file_name in listimg.ccds(ccd_kwargs={'unit':'adu'},return_fname=True):
        t120.log.info('now considering file '+file_name)
        listccd.append(ccd)
    
    combiner = Combiner(listccd)
    t120.log.info('now making the combination')
    master_flat = combiner.median_combine()
    fits_ccddata_writer(master_flat,master_file)
    t120.log.info('Result saved in '+master_file)
    return

def t120_mkflat(flat_dir=t120.t120_flat_dir,master_name_root='master',master_offset=t120.t120_ofst_dir+t120.t120_master_name):
    # read offset
    hdu_offset_list	= fits.open(master_offset)
    offset = CCDData(hdu_offset_list[0].data,unit=u.adu)
    
    dict_ccd_data = {}
    list_ccd_data = []
    for fit_file in glob.glob(flat_dir+'*.fit'):
     t120.log.info('now opening file: '+fit_file)
     hdu = fits.open(fit_file)
     filter_name = hdu[0].header['FILTER']
     t120.log.info('filter='+filter_name)
     if not dict_ccd_data.has_key(filter_name):
        dict_ccd_data[filter_name]=[]
     else:
        dict_ccd_data[filter_name].append(subtract_overscan(CCDData(hdu[0].data,unit=u.adu),offset))
    
    t120.log.info('now loop over the filters')
    for filter_name in dict_ccd_data:
        combiner = Combiner(dict_ccd_data[filter_name])
        master_flat = combiner.median_combine()
        hdu = master_flat.to_hdu()
        master_file = flat_dir+master_name_root+'-'+filter_name+'.fits'
        hdu.writeto(master_file,overwrite=True)
        t120.log.info('Master flat saved in '+master_file)
    return



def t120_makecosmetic(work_dir=t120.t120_data_path,subdir_root='/2017/2017-08-21/',
                      orig_dir_root=t120.t120_orig_dir,reduc_dir_root=t120.t120_redu_dir,
                      offset_dir_root=t120.t120_ofst_dir,flat_dir_root=t120.t120_flat_dir,
                      common_dir_root=t120.t120_common_root,
                      scampahead_file_name=t120.t120_scamp_ahead_root,
                      master_offset_name=t120.t120_master_name):
    orig_dir        = work_dir+subdir_root+orig_dir_root
    reduc_dir       = work_dir+subdir_root+reduc_dir_root
    offset_dir      = work_dir+subdir_root+offset_dir_root
    flat_dir        = work_dir+subdir_root+flat_dir_root
    common_dir      = work_dir+common_dir_root
    scampahead_file = common_dir+scampahead_file_name
    
    # make master offset
    t120.log.info('now making master offset')
    t120_mkoffset(offset_dir)
    # make master flat
    master_offset_file = offset_dir+master_offset_name
    t120.log.info('now making master flat')
    t120_mkflat(flat_dir,master_offset=master_offset_file)
    
    # now check existence of useful files and directories
    for subdir in [work_dir,work_dir+subdir_root,orig_dir,offset_dir,flat_dir,common_dir]:
        if not os.path.isdir(subdir):
            msg = '*** FATAL ERROR: directory '+subdir+' does not exist'
            t120.log.error(msg)
            raise IOError(msg)
    try:
        master_offset = fits_ccddata_reader(master_offset_file)
    except:
        msg = '*** FATAL ERROR while reading '+master_offset_file
        t120.log.error(msg)
        raise IOError(msg)
    
    # read scamp ahead file
    scamp_header = read_scamp_ahead(scampahead_file)
    
    listimg = ImageFileCollection(orig_dir+'/')
    for ccd,fit_file in listimg.ccds(ccd_kwargs={'unit':'adu'},return_fname=True,save_location=reduc_dir): #,save_with_name='-c',save_location=reduc_dir+'/'):
        t120.log.info('now treating file: '+orig_dir+fit_file)
        filter_name = ccd.header['FILTER']
        flat_name   = flat_dir+'/master-'+filter_name+'.fits'
        master_flat = fits_ccddata_reader(flat_name,unit=u.adu)
        master_offset.header['exposure'] = ccd.header['exposure'] # for dark subtraction
        ccd_corr    = ccd_process(ccd,exposure_key='exposure',
                                      exposure_unit=u.second,
                                      dark_frame=master_offset,
                                      master_flat=master_flat)
        ccd_corr.header             = put_scamp_header(ccd_corr.header,scamp_header)	# update header
        skycoo                      = SkyCoord(ccd_corr.header['OBJCTRA']+' '+ccd_corr.header['OBJCTDEC'],
                                        unit=(u.hourangle, u.deg))
        ccd_corr.header['RADESYS']  = ('ICRS','Astrometric system')
        ccd_corr.header['CUNIT1']   = ('deg','Axis unit')
        ccd_corr.header['CUNIT2']   = ('deg','Axis unit')
        ccd_corr.header['BUNIT']    = ('adu')
        ccd_corr.header['CRVAL1']   = (skycoo.ra.to('deg').value ,'Reference Right ascencion in decimal deg')
        ccd_corr.header['CRVAL2']   = (skycoo.dec.to('deg').value,'Reference Declination in decimal deg')
        ccd_corr.header['CTYPE1']   = ('RA---TAN','first parameter RA  ,  projection TANgential')
        ccd_corr.header['CTYPE2']   = ('DEC--TAN','second parameter DEC,  projection TANgential')
        out_fit_file                = reduc_dir+'/'+os.path.splitext(fit_file)[0]+'-c.fits'
        copy_fit_file               = reduc_dir+fit_file
        silentremove.silentremove(out_fit_file)
        silentremove.silentremove(copy_fit_file) # if exists produces error
        # make primary HDU
        hducorrlist			= ccd_corr.to_hdu()
        ccd_tosave			= CCDData(hducorrlist[0].data,unit=u.adu)
        ccd_tosave.header		= hducorrlist[0].header
        fits_ccddata_writer(ccd_tosave,out_fit_file)
        t120.log.info('corrected image saved in '+out_fit_file)
        #silentremove.silentremove(copy_fit_file)		# because it's useless...

"""
t120_makecosmetic(subdir_root='/2017/2017-08-21/')
t120_makecosmetic(subdir_root='/2017/2017-08-22/')
t120_makecosmetic(subdir_root='/2017/2017-08-23/')
"""
t120_makecosmetic(subdir_root='/2018/2018-08-13/')
