""" Launches astromatic programs for OHP T120 images.

See also
--------
www.astromatic.net
https://readthedocs.org/projects/sextractor/
https://readthedocs.org/projects/scamp/

"""

import os
import shutil
import datetime
import glob
import subprocess
import numpy as np
import astropy.units as u
from astropy.time import Time,TimeDelta
from astropy.io.votable import parse
from astropy.io import fits

from t120_init import t120


def put_all_new_header(path,mincontrast=2.0,astrodir=t120.t120_astr_dir):
    """Put SCAMP computed header into image files.
    
    Parameters
    ----------
    path : string
        root of the directory name where data are located.
    mincontrast : float, optional
        minimum contrast below which fits files are not updated.
    astrodir : string, optional
        Directory name where correct contrast files are copied.
    
    Returns
    -------
    None.
    
    """
    
    os.chdir(path)
    scampxml = path+t120.t120_scampres_dir + 'scamp.xml'
    if not os.path.isfile(scampxml):
        msg = '*** FATAL ERROR: file '+scampxml+' does not exist.'
        t120.log.error(msg)
        raise IOError(msg)
    # read XML file
    votable = parse(scampxml)
    table = votable.get_first_table()
    contrast = (table.array['XY_Contrast']).data
    catname = (table.array['Catalog_Name']).data
    
    for (contrast,catname) in zip((table.array['XY_Contrast']).data,(table.array['Catalog_Name']).data):
        if (contrast < mincontrast):
            img_file = path+t120.t120_redu_dir + catname.replace('.ldac','.fits')
            if not os.path.exists(img_file):
                t120.log.warning('*** WARNING: file '+img_file+' does not exists')
                continue
            header_file = path+t120.t120_scampres_dir + catname.replace('.ldac','.head')
            put_new_header(img_file,header_file)
            shutil.copyfile(img_file, img_file.replace(t120.t120_redu_dir,t120.t120_astr_dir))
    return
    
    

def put_new_header(img_file,header_file):
    """Put header in an image from a header file.
    
    The goal is to put astrometric information in a file that does not have it.
    Usually the header was computed by SCAMP.
    
    Parameters
    ----------
    img_file : string
        Fits file to update.
    header_file : string
        header file name containing the astrometry info. Usually it is computed by SCAMP.
        
    Returns
    -------
    None.
        The img_file is updated with the astrometry information contained in the head_file.
    
    
    """
    target = img_file
    # read image
    t120.log.info('reading image: '+img_file)
    hdul = fits.open(img_file)
    data_img = hdul[0].data
    hdr_img = hdul[0].header
    hdul.close()
    # clean image header
    #hdr_img = CleanHeaderNonAscii(hdr_img)
    # read header file
    t120.log.info('reading header: '+header_file)
    hdr_fil = fits.Header.fromfile(header_file,sep='\n',endcard=False,padding=False)
    # fix the infinite FLXSCALE possible problem
    #hdr_fil = fix_FLXSCALE(hdr_fil)
    # now loop over hdr_fil to updte the hdr_img
    for hdr in hdr_fil:
        hdr_img.set(hdr,value=hdr_fil.cards[hdr].value,comment=hdr_fil.cards[hdr].comment) #,verify('fix'))
    hdr_img.add_history('ASTROMETRY WAS TAKEN FROM FILE '+os.path.basename(header_file))
    # overwrite image
    fits.writeto(target,data_img,header=hdr_img,output_verify='exception', overwrite=True)
    t120.log.info(img_file+' was updated with astrometry info contained in '+header_file)
    return
    
def movepattern(patterns,outpath,strict=True):
    """Move the files in the pattern into the outpath.
    
    Parameters
    ----------
    pattern : array of string.
        Pattern of the files to move.
    outpath : string
        directory to put the files into.
    strict : boolean, optional
        If True, raise IOError if it is impossible to move the considered files.
    
    Returns
    -------
    None.
    
    """
    for pattern in patterns:
        listfile = glob.glob(pattern)
        t120.log.info('There are '+str(len(listfile))+' files of pattern= '+pattern)
        if len(listfile):
            for origf in listfile:
                destf = outpath+os.path.basename(origf)
                t120.log.debug('Now moving file '+origf+' into '+destf)
                try:
                    shutil.move(origf,destf)
                except:
                    msg = '*** FATAL ERROR: impossible to move '+origf+' into '+destf
                    t120.log.error(msg)
                    raise IOError(msg)
        else:
            t120.log.warning('There is no '+pattern+' file')
            if strict:
                raise IOError('There is no '+pattern+' file')
    return


def cmd_launch(cmd):
    """Submit a shell command"
    
    Raises exception if submission failed.
    
    Parameters
    ----------
    cmd : string
        shell command to submit
    
    Returns
    -------
    None.
    """
    t120.log.info('Launching cmd= '+''.join(cmd))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    t120.log.debug('out= ')
    t120.log.debug(out)
    t120.log.debug('err= ')
    t120.log.debug(err)
    return

def launch_astromatic(root,path='./',procdir=t120.t120_redu_dir,Lsex=False,Lpsfex=False,Lsexpsf=False,Lscamp=False,LscampNFOV=False,Lnewahead=False):
    """Launches SExtractor, PSFEX, SCAMP and update the .ahead file.
    
    Parameters
    ----------
    root : string
        root name of the input file.
    path : string, optional
        Path from where the script is launched, knowing that
        images are taken from the procdir directory. Default is './'.
    procdir : string, optional
        Processed directory containing corrected images to process.
        Default is t120.t120_redu_dir.
    Lsex : boolean, optional
        If True, SExtractor is launched.
    Lpsfex : boolean, optional
        If True, psfex is launched.
    Lsexpsf : boolean, optional
        If True, SExtractor is launched after psfex.
    Lscamp : boolean, optional
        If True, SCAMP is launched.
    Lnewahead : boolean, optional
        If True, the new ahead file is updated.
    LscampNFOV : boolean, optional
        If True, SCAMP is run twice: the 2nd time with a narrow field search for matching stars.
        Advantage: when there are many stars and SCAMP works great, the astrometry solution is better.
        The risk is that SCAMP cannot converge and make the solution worse.
        Default is False.
    
    Returns
    -------
    None.
    
    """
    rootsex='-c.fits'
    rootpsfex='-c.ldac'
    rootscamp='-c.ldac'
    refldac=''
    # first relocate to path
    Lchdir = False
    if not (path=='./'):
        Lchdir = True
    if Lchdir:
        curdir = os.getcwd()
        t120.log.debug('Relocating to '+path)
        os.chdir(path)
    #########################################################################
    # launches sextractor
    #########################################################################
    pattern = './'+t120.t120_redu_dir+root+'*'+rootsex
    listsex = glob.glob(pattern)
    if (Lsex):
        t120.log.info('=================================================')
        t120.log.info('============ Now lauching sextractor ============')
        t120.log.info('=================================================')
        if not len(listsex):
            msg = '*** FATAL ERROR: there is no file of type: '+pattern
            t120.log.warning(msg)
            raise IOError(msg)
        else:
            t120.log.info(' listsex='+' '.join(listsex))
            for filesex in listsex:
                t120.log.info(' --- Now treating file: '+filesex)
                filename = os.path.basename(filesex) # file without full path
                filebasename = os.path.splitext(filename)[0]  # file without extension
                catalogname = filesex.replace(rootsex,rootpsfex) # "filebasename.ldac"
                bckgname = filesex.replace(rootsex,'-bckg.fits')
                objname = filesex.replace(rootsex,'-obj.fits')
                cmd = [t120.t120_sex_exe+' '+filesex+' -c '+t120.t120_sex_conf+
                    ' -PARAMETERS_NAME '+t120.t120_sex_param+
                    ' -FILTER_NAME '+t120.t120_sex_filter+
                    #' -WEIGHT_IMAGE '+t120.t120_mask+
                    ' -CATALOG_NAME '+catalogname+
                    ' -CHECKIMAGE_NAME '+bckgname+','+objname]
                cmd_launch(cmd[0])
                # now moves the sextractor results into $t120.t120_sexres_dir subdirectory
                sexresdir = './'+t120.t120_sexres_dir
                if not os.path.exists(sexresdir):
                    t120.log.debug('Now making directory '+sexresdir)
                    os.mkdir(sexresdir)
                #movepattern([catalogname,bckgname,objname],sexresdir)
                movepattern([catalogname],sexresdir)
    #########################################################################
    # launches PSFEX
    #########################################################################
    if (Lpsfex):
        t120.log.info('=================================================')
        t120.log.info('============ Now lauching psfex ============')
        t120.log.info('=================================================')
        pattern = './'+t120.t120_sexres_dir+root+'*'+rootpsfex
        listpsfex = glob.glob(pattern)
        if not len(listpsfex):
            msg = '*** FATAL ERROR: there is no file of type: '+pattern
            t120.log.warning(msg)
            raise IOError(msg)
        else:
            cmd = t120.t120_psfex_exe+' -c '+t120.t120_psfex_conf+' '+' '.join(listpsfex)
            cmd_launch(cmd)
            # now moves psfex results into PSF dir
            psfexdir = './'+t120.t120_psfres
            if not os.path.exists(psfexdir):
                t120.log.debug('Now making directory '+psfexdir)
                os.mkdir(psfexdir)
            patterns = ['chi*','proto*','samp*','snap*','resi*','psfex.xml',
                    'countfrac*','counts_*','ellipticity_*','fwhm_*']
            movepattern(patterns,psfexdir)
    #########################################################################
    # launches SExtractor after PSFEX
    #########################################################################
    if (Lsexpsf):
        t120.log.info(' =============================================================')
        t120.log.info(' ============ Now RE-lauching sextractor with psf ============')
        t120.log.info(' =============================================================')
        pattern = './'+t120.t120_redu_dir+root+'*'+rootsex
        listsex = glob.glob(pattern)
        if not len(listsex):
            msg = '*** FATAL ERROR: there is no file of type: '+pattern
            t120.log.warning(msg)
            raise IOError(msg)
        else:
            t120.log.info(' listsex='+' '.join(listsex))
            for filesex in listsex:
                t120.log.info(' --- Now treating file: '+filesex)
                filename = os.path.basename(filesex) # file without full path
                filebasename = os.path.splitext(filename)[0]  # file without extension
                psfile = t120.t120_psfres + filesex.replace(rootsex,'psf')
                catalogname = filesex.replace(rootsex,rootpsfex) # "filebasename.ldac"
                bckgname = filesex.replace(rootsex,'-bckg.fits')
                objname = filesex.replace(rootsex,'-obj.fits')
                cmd = [t120.t120_sex_exe+' '+filesex+' -c '+t120.t120_sex_conf+
                    ' -PSF_NAME '+psfile+
                    ' -PARAMETERS_NAME '+t120.t120_sex_parampsf+
                    ' -FILTER_NAME '+t120.t120_sex_filter+
                    ' -WEIGHT_IMAGE '+t120.t120_mask+
                    ' -CATALOG_NAME '+catalogname+
                    ' -CHECKIMAGE_NAME '+bckgname+','+objname]
                cmd_launch(cmd[0])
                # now moves the sextractor results into $t120.t120_sexres_dir subdirectory
                sexresdir = './'+t120.t120_sexres_dir
                if not os.path.exists(sexresdir):
                    t120.log.debug('Now making directory '+sexresdir)
                    os.mkdir(sexresdir)
                movepattern(['*.ldac','bckgname','objname'],sexresdir)
    #########################################################################
    # launches SCAMP
    #########################################################################
    for (llscamp,scamp_conf_file) in zip([Lscamp,LscampNFOV],[t120.t120_scamp_conf_WFOV,t120.t120_scamp_conf_NFOV]):
        if (llscamp):
            pattern = t120.t120_sexres_dir+root+'*'+rootscamp
            t120.log.debug('pattern: '+pattern)
            listscamp = glob.glob(pattern)
            if not len(listscamp):
                msg = '*** FATAL ERROR: there is no file of type: '+pattern
                t120.log.warning(msg)
                raise IOError(msg)
            else:
                t120.log.debug('listscamp='+' '.join(listscamp))
                t120.log.info('============================================')
                t120.log.info('============ Now lauching scamp ============')
                t120.log.info('============================================')
                cmd = [t120.t120_scamp_exe+' '+ ' '.join(listscamp)+
                    ' -c '+scamp_conf_file] #+' -AHEADER_GLOBAL '+t120.t120_scamp_ahead]
                cmd_launch(cmd[0])
                # now move SCAMP results into scamp dir
                scampresdir = './'+t120.t120_scampres_dir
                if not os.path.exists(scampresdir):
                    t120.log.debug('Now making directory '+scampresdir)
                    os.mkdir(scampresdir)
                movepattern(['*.png','scamp.xml',t120.t120_sexres_dir+'*.head'],scampresdir)
                # now moves the correct astrometry images into astrometry directory
                # first open the scamp.xmk file
    """
    #########################################################################
    # Determine best contrast file and copy it to default
    #########################################################################
    if (Lnewahead):
        t120.log.info('==================================================================')
        t120.log.info('============ Now determining best contrast file ==================')
        t120.log.info('==================================================================')
        bestcontfile = get_bestcontrast()
        t120.log.info('bestcontfile='+bestcontfile)
        head2ahead(bestheadfile,camera.ahead_file)
    """
    #########################################################################
    #########################################################################
    # go back to initial directory
    if Lchdir:
        t120.log.debug('Relocating to '+curdir)
        os.chdir(curdir)
    
    # end of program
    t120.log.info('done')
    return

