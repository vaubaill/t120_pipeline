import glob
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
import t120_getradec

path = '/Volumes/VAUBGO/2018-08-13-DU-ECU/REDUCED/'
root = '21P'
loc='511'

pattern = path+root+'*.fits'
fitsfile=glob.glob(pattern)[0]
(objname,time) = t120_getradec.get_nametime(fitsfile)
print(fitsfile+' '+objname+' '+time.isot)
skyc = t120_getradec.get_radec_fromjpl(ssid=objname,location='511',jd=time.jd)

