# T120 pipeline initialization
#
# J. Vaubaillon, IMCCE, 2018
#
import os
import sys
import glob
import shutil
import logging
import subprocess
import numpy as np
from configparser import ConfigParser, ExtendedInterpolation
import astropy.units as u
from astropy.io.votable import parse
from astropy.io import fits
from astropy.time import Time as Time
from astropy.coordinates import SkyCoord, Angle, Latitude, Longitude, EarthLocation, AltAz, ICRS
from astropy.table import QTable
import astropy.io.ascii as ascii

from t120_log import get_module_logger,get_log_file_name

home=os.environ['HOME']+'/'

class T120(object):
    """T120 object definition.
    
    """
    config_file='../common/t120.ini'
    config = ConfigParser(interpolation=ExtendedInterpolation())
    config.read(config_file)
    # home directory
    t120_home        = home
    # T120Pipeline full path
    t120_pipe_path   = config['PATH']['pipe_path']
    # T120 data (sub)-directory
    t120_data_path = config['PATH']['data_path']
    # T120 common dir root
    t120_common_root = config['PATH']['common_root'] 
    # T120 common directory where all external configuration files are located
    t120_common_dir = config['PATH']['common_dir'] 
    # T120 ORIGINAL data directory name
    t120_orig_dir = config['PATH']['orig_dir']
    # T120 REDUCED data directory name
    t120_redu_dir = config['PATH']['redu_dir']
    # T120 FLAT data directory name
    t120_flat_dir = config['PATH']['flat_dir']
    # T120 DARK data directory name
    t120_dark_dir = config['PATH']['dark_dir']
    # T120 OFFSET data directory name
    t120_ofst_dir = config['PATH']['ofst_dir']
    # T120 master offst name
    t120_master_name = config['PATH']['master_name']
    # T120 Astrometry data directory name
    t120_astr_dir = config['PATH']['astr_dir']
    
    
    # T120 JPL id text file name
    t120_idjplfile = config['JPL']['idjplfile']
    # T120 id radec text file name
    t120_idradecfile = config['JPL']['idradecfile']
    
    # Astromatic
    # SExtractor executable
    t120_sex_exe = config['ASTROMATIC']['sex_exe']
    # SExtractor configuration file
    t120_sex_conf = config['ASTROMATIC']['sex_conf']
    # SExtractor parameters file
    t120_sex_param = config['ASTROMATIC']['sex_param']
    # SExtractor parameters given PSFEX file
    t120_sex_parampsf = config['ASTROMATIC']['sex_parampsf']
    # SExtractor filter file
    t120_sex_filter = config['ASTROMATIC']['sex_filter']
    # SCAMP executable
    t120_scamp_exe = config['ASTROMATIC']['scamp_exe']
    # SCAMP configuration file for wide FOV
    t120_scamp_conf_WFOV = config['ASTROMATIC']['scamp_conf_WFOV']
    # SCAMP configuration file for narrow FOV
    t120_scamp_conf_NFOV = config['ASTROMATIC']['scamp_conf_NFOV']
    # PSFEX executable file
    t120_psfex_exe = config['ASTROMATIC']['psfex_exe']
    # PSFEX configuration file
    t120_psfex_conf = config['ASTROMATIC']['psfex_conf']
    # MISSFITS executable
    t120_missfits_exe = config['ASTROMATIC']['missfits_exe']
    # MISSFITS configuration file
    t120_missfitsconfig = config['ASTROMATIC']['missfitsconfig']
    
    # T120 scamp ahead file name
    t120_scamp_ahead_root = config['ASTROMETRY']['scamp_ahead_root']
    # T120 scamp ahead file full path
    t120_scamp_ahead = config['ASTROMETRY']['scamp_ahead']
    # T120 SCAMP Resutls directory name
    t120_scampres_dir = config['ASTROMETRY']['scampres_dir']
    # T120 SExtractor Resutls directory name
    t120_sexres_dir = config['ASTROMETRY']['sexres_dir']
    
    # logging directory
    t120_log_path = config['PATH']['log_path']
    # default log file
    log_dft_file = get_log_file_name(pathlog=t120_log_path)
    log = get_module_logger('t120',logfile=log_dft_file,addsldr=True)
    
    t120_toprocdir = config['TOPROCESS']['procdir']

# creates a CAB object
t120=T120()




def makeobslog(path,root=''):
    """Make a journal of an observation night.
    
    Parameters
    ----------
    path : string
        path where data are located.
    root : string, optional
        root of files to journal.
        Default is '', i.e. all files will be journaled.
    
    Returns
    -------
    None.
    
    """
    
    # user definitions example
    # path = t120.t120_data_path + '2019/2019-07-29-DU-ECU/ORIGINAL/'
    # root = ''
    
    # start of program: do not touch
    fileout = path+ 'journal' + root + '.log'
    pattern = path +  root + '*.fit*'
    listfile = glob.glob(pattern)
    nfile = len(listfile)
    if (nfile==0):
        msg = '*** WARNING: there is no file of type: '+pattern
        t120.log.warning(msg)
        raise IOError(msg)
    t120.log.info('There are '+str(nfile)+' files of type '+pattern)
    prefix = os.path.commonprefix(listfile)
    # create output table: make a first dummy row so astropy.table.QTable knows 'FILTER' is a string
    data = QTable([['                                            '],['                    '],
                [0.0],[0.0],[0.0],[0.0],[' '],[0.0]],
                names=['FILE','TARGET','JD','RA (deg)','DEC (deg)','EXP (s)','FILTER','AIRMASS'])
    # loop over the files
    for imgfile in listfile:
        hdul = fits.open(imgfile)
        hdr = hdul[0].header
        target_name = hdr['OBJECT']#.strip().upper().replace(' ','')
        # JD, the center coordinate (RA, Dec), exposure time, filter, airmass.
        JD = hdr['JD']
        try:
            RAcenter = hdr['CRVAL1']
            DEcenter = hdr['CRVAL2']
        except:
            t120.log.warning('There is no CRVAL keyword  in file '+imgfile)
            skycoo = SkyCoord(hdr['OBJCTRA'],hdr['OBJCTDEC'], unit=(u.hourangle, u.deg))
            RAcenter = skycoo.ra.to('deg').value
            DEcenter = skycoo.dec.to('deg').value
        exptime = hdr['EXPOSURE']
        filtr = hdr['FILTER']
        airmass = hdr['AIRMASS']
        data.add_row([imgfile.replace(prefix,'').lstrip(),target_name,JD,RAcenter,DEcenter,exptime,filtr,airmass])
        t120.log.debug('File='+imgfile+' TARGET'+target_name+' JD='+str(JD)+' RA='+str(RAcenter)+' DE='+str(DEcenter)+
                  ' exp='+str(exptime)+' fil='+str(filtr)+' airmass='+str(airmass))
    # remove first row
    data.remove_row(0)
    # save log in file
    ascii.write(data,fileout,format='fixed_width',delimiter=' ',formats={'JD': '%18.12f'},overwrite=True)
    t120.log.info('There are '+str(len(listfile))+' data saved in '+fileout)
    return
