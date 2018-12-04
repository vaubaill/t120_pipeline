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
import astropy.units as u
from astropy.io.votable import parse
from astropy.io import fits
from astropy.time import Time as Time
from astropy.coordinates import SkyCoord, Angle, Latitude, Longitude, EarthLocation, AltAz, ICRS

home=os.environ['HOME']+'/'

class T120(object):
    """T120 object definition.
    
    """
    # home directory
    t120_home        = home
    # T120Pipeline full path
    t120_pipe_path   = home + '/OBSERVATIONS/OHPT120/python/'
    # T120 data (sub)-directory
    #t120_data_path   = '/Volume/BackUpMac/OHP/'
    t120_data_path = '/Users/vaubaill/OBSERVATIONS/OHPT120/2018/'
    # T120 common dir root
    t120_common_root = 'common/'
    # T120 common directory where all external configuration files are located
    t120_common_dir = t120_data_path + t120_common_root
    # T120 ORIGINAL data directory name
    t120_orig_dir = 'ORIGINAL/'
    # T120 REDUCED data directory name
    t120_redu_dir = 'REDUCED/'
    # T120 FLAT data directory name
    t120_flat_dir = 'FLAT/'
    # T120 DARK data directory name
    t120_dark_dir = 'DARK/'
    # T120 OFFSET data directory name
    t120_ofst_dir = 'OFFSET/'
    # T120 master offst name
    t120_master_name = 'master.fits'
    # T120 Astrometry data directory name
    t120_astr_dir = 'Astrometry/'
    
    
    # T120 JPL id text file name
    t120_idjplfile = 'nameidjpl.txt'
    # T120 id radec text file name
    t120_idradecfile = 'idradec.txt'
    
    # T120 scamp ahead file name
    t120_scamp_ahead_root = 't120-andor1024.ahead'
    # T120 scamp ahead file full path
    t120_scamp_ahead = t120_common_dir + t120_scamp_ahead_root
    # T120 SCAMP Resutls directory name
    t120_scampres_dir = 'SCAMPRES/'
    
    # T120 SExtractor Resutls directory name
    t120_sexres_dir = 'SEXRES/'
    
    # logging directory
    t120_log_path    = t120_pipe_path + 'LOG/'
    # default log file
    log_dft_file = t120_log_path + 't120-' + Time.now().isot.replace('-','').replace(':','').split('.')[0] + '.log' 
    #log_dft_file = t120_log_path + 'CABPIPE.log' 
    # log format
    log_fmt = '%(levelname)s %(filename)s %(lineno)d (%(funcName)s) : %(message)s '
    # logger name
    log_name='cab'
    # create logger
    log = logging.getLogger(log_name)
    # set logger level
    log.setLevel(logging.DEBUG)
    log.propagate = True
    # create StreamHandler
    sdlr = logging.StreamHandler()
    # set StreamHandler options
    sdlr.setLevel(logging.INFO)
    sdlr.setFormatter(logging.Formatter(fmt=log_fmt))
    # add StreamHandler into logger
    log.addHandler(sdlr)
    # creates a default file handler
    hdlr = logging.FileHandler(log_dft_file)
    # set file handler options
    hdlr.setLevel(logging.DEBUG)
    hdlr.setFormatter(logging.Formatter(fmt=log_fmt))
    # adds file handler to logger
    log.addHandler(hdlr)
    log.info('Default log file is: '+log_dft_file)

# creates ca CAB object
t120=T120()

