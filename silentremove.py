import os, errno

def silentremove(filename):
    try:
        os.remove(filename)
	print 'file ',filename,' has been removed'
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        print 'unable to remove ',filename
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred
