import glob
from astropy.io import fits
from astropy.time import Time
from t120_init import t120
import t120_getradec

#path = '/Volumes/VAUBGO/2018-08-13-DU-ECU/REDUCED/'
path = '/Users/vaubaill/OBSERVATIONS/OHPT120/2018/2018-08-15-DU-ECU/'

root = 'REDUCED/P10J8BV'
loc='511'

Ljpl=True
Lname=True
Ljpltxt=True
Lnametxt=True

blockJPL=False
blockName=False
blockJPLtxt=False
blockTxt=False
blockFound=True

pattern = path+root+'*.fits'
listfile=glob.glob(pattern)
if (len(listfile)==0):
    raise IOError('*** FATAL ERROR: no file of type '+pattern)
for fitsfile in listfile:
    skycoo = None
    (objname,time) = t120_getradec.get_nametime(fitsfile)
    t120.log.info(fitsfile+' '+objname+' '+time.isot)
    # try to get RA,DEC from astroquery.horizons
    if Ljpl:
        t120.log.info('\n now trying to get RADEC of '+objname+' from astroquery.horizons')
        try:
            skycoo = t120_getradec.get_radec_fromjpl(ssid=objname,location=loc,jd=time.jd)
        except Exception:
            msg = '*** WARNING: impossible to get RA,DEC from JPL.'
            t120.log.warning(msg)
            if blockJPL:
                raise ValueError(msg)
    # try to get RA,DEC from astropy.coordinate.SkyCoord().from_name() function
    if Lname:
        if not skycoo:
            t120.log.info('\n now trying to get RADEC of '+objname+' with astropy.SkyCoord')
            try:
                skycoo = t120_getradec.get_radec_fromname(objname)
            except Exception:
                msg = '*** WARNING: impossible to get RA,DEC from name: '+objname
                t120.log.warning(msg)
                if blockName:
                 raise ValueError(msg)
    # try to get RA,DEC from a text file providing JPL id from Object name
    if Ljpltxt:
        if not skycoo:
            txtfile = path + t120.t120_idjplfile
            t120.log.info('\n now trying to get RADEC of '+objname+' from file '+txtfile)
            try:
                ssjplid = t120_getradec.get_jplid_fromtxtid(objname,txtfile)
                t120.log.info('ssjplid='+ssjplid)
                skycoo = t120_getradec.get_radec_fromjpl(ssid=ssjplid,location=loc,jd=time.jd)
            except Exception:
                msg = '*** WARNING: impossible to get RA,DEC from JPL id file: '+txtfile
                t120.log.warning(msg)
                if blockJPLtxt:
                    raise ValueError(msg)
    # try to get RA,DEC from a text file RA,DEC from Object name
    if Lnametxt:
        if not skycoo:
            txtfile = path + t120.t120_idradecfile
            t120.log.info('\n now trying to get RADEC of '+objname+' from file '+txtfile)
            try:
                skycoo = t120_getradec.get_radec_fromtxt(objname,txtfile)
            except Exception:
                msg = '*** WARNING: impossible to get RA,DEC from id RADEC file: '+txtfile
                t120.log.warning(msg)
                if blockTxt:
                    raise ValueError(msg)
    if not skycoo:
        msg = '\n \n *** FATAL ERROR: Impossible to get automated RA,DEC for object: '+objname
        t120.log.error(msg)
        t120.log.error('File '+fitsfile+' was NOT updated with RA,DEC.')
        if blockFound:
            raise ValueError(msg)
    else:
        t120_getradec.put_radec(fitsfile,skycoo)

t120.log.info('Done.')
