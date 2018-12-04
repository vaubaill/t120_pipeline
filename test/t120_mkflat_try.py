import glob
from astropy import units as u
from astropy.io import fits
from ccdproc import CCDData, Combiner
import numpy as np
from ccdproc import Combiner, ImageFileCollection,fits_ccddata_writer
from t120_init import t120

listimg = ImageFileCollection(t120.t120_flat_dir) #,glob_include='*.fit',glob_exclude='*.fits')

list_filters = listimg.values('filter',unique=True)
#for filter_name in listimg.values('filter',unique=True):
for filter_name in list_filters:
    t120.log.info('*** filter: '+filter_name)
    my_files = listimg.files_filtered(filter=filter_name)
    t120.log.info('my_files='+my_files)
 
pouet


for filter_name in listimg.values('filter',unique=True):
    t120.log.info('*** filter: '+filter_name)
    listccd = []
    for ccd,file_name in listimg.ccds(ccd_kwargs={'unit':'adu'},filter=filter_name,return_fname=True):
        t120.log.info('now considering file '+file_name)
        listccd.append(ccd)
    """
    t120.log.info('now making the Combiner object')
    combiner = Combiner(listccd)
    t120.log.info('now making the combination')
    master_flat = combiner.average_combine()
    master_file = t120.t120_flat_dir+'/master-'+filter_name+'.fits'
    fits_ccddata_writer(master_flat,master_file)
    t120.log.info('Master flat saved in '+master_file)
    """
