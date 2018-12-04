#import t120_makecosmetic
from t120_init import t120
import t120_getradec
import t120_cleanastrodir

# make journal of the night
#t120_mksummary(journal)
# make offset
#t120_mkoffset(journal=journal,path='OFFSET/')
# make dark
#t120_mkdark,journal=journal)
# make flat
#t120_mkflat,journal=journal)
# performs cosmetic operations
#t120_mkcosmetic,scampaheadfile=t120.t120_scamp_ahead_root)

# get RA,DEC of FOV
#t120_getradec.t120_getradec()

# launches SExtractor and SCAMP with large values of POSANGLE_MAXERR = 5.0 and POSITION_MAXERR = 15.0
#t120_sexscamp()
# relaunches SExtractor and SCAMP with small values of POSANGLE_MAXERR = 1.0 and POSITION_MAXERR = 5.0
#t120_sexscamp()
# clean Astrometry directory
t120_cleanastrodir()

#spawn,"pushd ./ ; cd /Users/vaubaill/OBSERVATIONS/OHPT120/python ; python test_getradec_fromjpl.py ; popd"
#spawn,"sext120.sh"					# launches Astrometry suite

#t120_cleanallhdr
