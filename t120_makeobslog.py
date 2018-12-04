"""Make observation log.



"""
import os
import glob
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.table import QTable
import astropy.io.ascii as ascii

from t120_init import t120

# user definitions
path = t120.t120_data_path + '2018*/Astrometry/'
root = 'Phaethon'

# start of program: do not touch
fileout = t120.t120_data_path + root + '-log.txt'
pattern = path +  root + '*.fits'
listfile = glob.glob(pattern)
nfile = len(listfile)
if (nfile==0):
    msg = '*** WARNING: there is no file of type: '+pattern
    t120.log.warning(msg)
    raise IOError(msg)
t120.log.info('There are '+str(nfile)+' files of type '+pattern)
prefix = os.path.commonprefix(listfile)
# create output table: make a first dummy row so astropy.table.QTable knows 'FILTER' is a string
data = QTable([['                                            '],
                [0.0],[0.0],[0.0],[0.0],[' '],[0.0]],
                names=['FILE','JD','RA (deg)','DEC (deg)','EXP (s)','FILTER','AIRMASS'])
# loop over the files
for imgfile in listfile:
    hdul = fits.open(imgfile)
    hdr = hdul[0].header
    # JD, the center coordinate (RA, Dec), exposure time, filter, airmass.
    JD = hdr['JD']
    RAcenter = hdr['CRVAL1']
    DEcenter = hdr['CRVAL2']
    exptime = hdr['EXPOSURE']
    filtr = hdr['FILTER']
    airmass = hdr['AIRMASS']
    data.add_row(['201'+imgfile.replace(prefix,'').lstrip(),JD,RAcenter,DEcenter,exptime,filtr,airmass])
    t120.log.debug('File='+imgfile+' JD='+str(JD)+' RA='+str(RAcenter)+' DE='+str(DEcenter)+
                  ' exp='+str(exptime)+' fil='+str(filtr)+' airmass='+str(airmass))
# remove first row
data.remove_row(0)
# save log in file
ascii.write(data,fileout,format='fixed_width',delimiter=' ',formats={'JD': '%18.12f'},overwrite=True)
t120.log.info('There are '+str(len(listfile))+' data saved in '+fileout)
