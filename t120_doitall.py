import os
import logging
import glob

from t120_init import t120,makeobslog
from t120_clean import cleanastrodir,cleancosmeticdir
from t120_getradec import t120_getradec
from t120_insert_radec import t120_insert_radec
from t120_makecosmetic import t120_make_all_cosmetic,t120_makecosmetic,t120_mkoffset,t120_mkflat,t120_mkdark
from t120_astromatic import launch_astromatic, put_all_new_header


allsubpath = glob.glob(t120.t120_toprocdir+'*')
t120.log.info('allsubpath '+' '.join(allsubpath))
for subpath in allsubpath:
    subpath= subpath+'/'
    t120.log.info('Now treating directory: '+subpath)
    """
    # make journal
    makeobslog(subpath+'ORIGINAL/')
    # first clean directories
    cleancosmeticdir(subpath)
    # set names
    offset_dir  = subpath+t120.t120_ofst_dir
    flat_dir    = subpath+t120.t120_flat_dir
    # make master offset
    t120.log.info('now making master offset')
    master_offset_file=t120_mkoffset(offset_dir)
    # make master flat
    t120.log.info('now making master flat')
    t120_mkflat(flat_dir,master_offset=master_offset_file)
    offset_dir  = subpath+t120.t120_ofst_dir
    offset_file  = offset_dir+t120.t120_master_name
    dark_dir = subpath+t120.t120_dark_dir
    t120_mkdark(dark_dir=dark_dir,master_offset=offset_file)
    # make cosmetic corrections
    t120_makecosmetic(work_dir=subpath)
    # put approximate radec into files
    t120_insert_radec(subpath)
    # put more accurate RADEC into file
    t120_getradec(path=subpath+t120.t120_redu_dir)
    """
    # makes astrometry
    launch_astromatic('',path=subpath,procdir='REDUCED/',Lsex=False,Lscamp=True,LscampNFOV=False)
    # put new header into files
    put_all_new_header(subpath)
    
    """
    # clean Astrometry directory
    t120_cleanastrodir(subpath)
    """


#t120_sexscamp()
# relaunches SExtractor and SCAMP with small values of POSANGLE_MAXERR = 1.0 and POSITION_MAXERR = 5.0
#t120_sexscamp()

#spawn,"pushd ./ ; cd /Users/vaubaill/OBSERVATIONS/OHPT120/python ; python test_getradec_fromjpl.py ; popd"
#spawn,"sext120.sh"					# launches Astrometry suite

#t120_cleanallhdr
