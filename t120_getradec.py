"""Get object ephemerides.


"""
import os
import glob
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy import units as u
from t120_init import t120
from astroquery.jplhorizons import conf, Horizons
from astropy.io import fits

def get_name(fitsfile):
    """Extract Solar system object from fits file using header.
    
    Parameters
    ----------
    fitsfile : string
        fits file name.
    
    Returns
    -------
    ssolname : string
        Solar System object name.
    
    """
    hdul = fits.open(fitsfile)
    hdr = hdul[0].header
    ssolname = hdr['OBJECT'].upper().replace(' ','')
    return ssolname

def get_nametime(fitsfile):
    """Extract Solar system object name and time from fits file using header.
    
    Parameters
    ----------
    fitsfile : string
        fits file name.
    
    Returns
    -------
    ssolname : string
        Solar System object name.
    time : astropy.time.Time Object
        Time of the observation.
    
    """
    hdul = fits.open(fitsfile)
    hdr = hdul[0].header
    ssolname = hdr['OBJECT'].upper().strip()
    time = Time(hdr['DATE-OBS'],format='isot')
    hdul.close()
    return (ssolname,time)

def defhorserver():
    """Change astroquery.jplhorizons (HORIZONS) server config.
    
    Parameters
    ----------
    None.
    
    Returns
    -------
    None.
    
    """
    horserver = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'
    conf.horizons_server = horserver
    t120.log.debug('HORIZONS server changed in to '+horserver)
    return
    
def get_radec_fromjpl(ssid='',location='511',jd=2454545.0):
    """Get Solar System Object position from Horizons
    
    Parameters
    ----------
    ssid : string, optional
        JPL id of the object to be found.
        Default is empty string.
    location : string, optional
        IAU observatory code. Default is '511' (OHP).
    jd : float, optional
        Juliand day for the computation of the ephemeris.
        Default is J2000.
    
    Returns
    -------
    skycoo : astropy.coordinates.Skycoordinate Object
        Sky coordinates (RA,DEC) of the object.
    
    """
    defhorserver()
    sshor = Horizons(id=ssid,location=location,epochs=jd)
    try:
        eph = sshor.ephemerides()
    except Exception as e:
        if ('Ambiguous' in str(e)):
            ssidjpl = str(e).split('\n')[-2].split()[0]
            t120.log.warning('Ambiguous name detected! Corrected to id='+ssidjpl)
            try:
                sshor = Horizons(id=ssidjpl,location=location,epochs=jd)
                eph = sshor.ephemerides()
            except Exception as e:
                raise(e)
        else:
            raise(e)
    skycoo = SkyCoord(eph['RA'][0],eph['DEC'][0],unit='deg',obstime=Time(jd,format='jd'))
    t120.log.info(ssid+' '+skycoo.to_string(style='hmsdms'))
    return skycoo


def get_radec_fromname(name):
    """Get sky object coordinate from its name.
    
    Note: NOT to be used by Solar system Object.
    It uses the astropy.coordinates.SkyCoord.from_name() method.
    
    Parameters
    ----------
    name : string
        name of the object
    
    Returns
    -------
    skycoo : astropy.coordinates.SkyCoord Object
        Sky coordinates of the object.
    
    """
    skycoo = SkyCoord(0,0,unit='deg').from_name(name)
    t120.log.info(name+' found by astropy.coordinates.SkyCoord')
    return skycoo 

def get_radec_fromtxt(name,txt):
    """Get sky coordinate from a text file.
    
    If name is not found in txt, None is returned.
    Note: in the file, object coordinates MUST be in the format:
    '12:34:56 +12:34:56' or '12:34:56.78 +12:34:56.78'
    i.e. column mark separated, NOT blank separated!
    
    Parameters
    ----------
    name : string
        Object name.
    txt : string
        File name containing the list of name and RA,DEC.
    
    
    Returns
    -------
    skycoo : astropy.coordinates.SkyCoord Object
        Sky coordinates of the object.
    
    
    """
    # read txt file
    skycoo = None
    if not os.path.exists(txt):
        msg = '*** FATAL ERROR: file '+txt+' does not exist'
        t120.log.error(msg)
        raise IOError(msg)
    found = False
    with open(txt,'r') as f:
        for line in f.readlines():
            if name.upper() in line.upper():
                radec = line.split("'")[1]
                skycoo = SkyCoord(radec.split()[0],radec.split()[1], unit=(u.hourangle, u.deg))
                found = True
    if not found:
        msg = '*** FATAL ERROR: '+name+' was not found in '+txt
        t120.log.error(msg)
        raise IOError(msg)
    
    return skycoo 

def get_jplid_fromtxtid(name,txtid):
    """
    
    Parameters
    ----------
    name : string
        Object name.
    txtid : string
        File name containing the list of name and JPL ids.
    
    Returns
    -------
    jplid : string
        JPL id.
    
    """
    if not os.path.exists(txtid):
        msg = '*** FATAL ERROR: file '+txtid+' does not exist'
        t120.log.error(msg)
        raise IOError(msg)
    # read txtid file
    found = False
    t120.log.info('Now reading file: '+txtid+' and searching for '+name)
    with open(txtid,'r') as f:
        for line in f.readlines():
            if name.upper() in line.upper():
                t120.log.info(name.upper()+' found in file: '+txtid)
                jplid = line.split()[1]
                found = True
    if not found:
        msg = '*** FATAL ERROR: '+name+' was not found in '+txtid
        t120.log.error(msg)
        raise IOError(msg)
    return jplid

def put_radec(fitsfile,skycoo):
    """Put RA,DEC coordinates into a fits file.
    
    Parameters
    ----------
    fitsfile : string
        Fits file name to update.
    skycoo : astropy.coordinates.SkyCoord Object
        Sky coordinates of the object.
    """
    comm = {'OBJCTRA':'Nominal Right Ascension of center of image',
            'OBJCTDEC':'Nominal Declination of center of image'}
    
    hdulist = fits.open(fitsfile)
    hdu = hdulist[0]
    data = hdu.data
    try:
        target_name = hdu.header['OBJECT'].strip()
    except:
        msg = '*** FATAL ERROR: target '+target_name+' not found'
        t120.log.error(msg)
        raiseValueError(msg)
    hdu.header['OBJECT'] = (target_name.upper().replace(' ',''),'Target name')
    t120.log.info('target_name='+target_name)
    time = Time(hdu.header['DATE-OBS'])
    site = EarthLocation(lat=hdu.header['SITELAT'], lon=hdu.header['SITELONG']) 
    altazframe = AltAz(obstime=time, location=site)
    altaz = skycoo.transform_to(altazframe)
    hdu.header['OBJCTRA'] = (skycoo.ra.to('deg').value,comm['OBJCTRA'])
    hdu.header['OBJCTDEC'] = (skycoo.dec.to('deg').value,comm['OBJCTDEC'])
    hdu.header['CRVAL1'] = (skycoo.ra.to('deg').value ,'Reference Right ascencion in decimal deg')
    hdu.header['CRVAL2'] = (skycoo.dec.to('deg').value,'Reference Declination in decimal deg')
    hdu.header['HISTORY'] = 'RA,DEC values inserted with py routine.'
    hdu.writeto(fitsfile,overwrite=True)
    t120.log.info('File '+fitsfile+' successfully updated with RADEC info')
    return


def t120_getradec(path=t120.t120_data_path+t120.t120_redu_dir,
        root = '*',loc='511',
        Ljpl=True,Lname=True,Ljpltxt=True,Lnametxt=True,
        blockJPL=False,blockName=False,blockJPLtxt=False,blockTxt=False,blockFound=True):
    """Get RA,DEC of center of FOV from object coordinates.
    
    TBF
    
    Parameters
    ----------
    TBF
    
    
    
    Returns
    -------
    TBF
    
    """
    pattern = path+root+'*-c.fits'
    listfile=glob.glob(pattern)
    if (len(listfile)==0):
        raise IOError('*** FATAL ERROR: no file of type '+pattern)
    for fitsfile in listfile:
        if '-c-c.fits' in fitsfile:
            msg='*** Warning: file with double -c will be deleted: '+fitsfile
            os.remove(fitsfile)
            continue
        skycoo = None
        (objname,time) = get_nametime(fitsfile)
        t120.log.info('*** Now processing: '+fitsfile+' '+objname+' '+time.isot)
        # try to get RA,DEC from astroquery.horizons
        if Ljpl:
            t120.log.info('\n now trying to get RADEC of '+objname+' from astroquery.horizons')
            try:
                skycoo = get_radec_fromjpl(ssid=objname,location=loc,jd=time.jd)
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
                    skycoo = get_radec_fromname(objname)
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
                    ssjplid = get_jplid_fromtxtid(objname,txtfile)
                    t120.log.info('ssjplid='+ssjplid)
                    skycoo = get_radec_fromjpl(ssid=ssjplid,location=loc,jd=time.jd)
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
                    skycoo = get_radec_fromtxt(objname,txtfile)
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
            t120.log.info('\n \n *** '+objname+' coordinates at '+time.isot+' : '+skycoo.to_string())
            put_radec(fitsfile,skycoo)
        
    t120.log.info('Done.')
    return

