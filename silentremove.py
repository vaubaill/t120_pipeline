import os, errno

def silentremove(filename):
    try:
        os.remove(filename)
        print 'file ',filename,' has been removed'
    except IOError as a:
        msg = '*** WARNING: unable to remove '+filename
        raise IOError(msg)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        msg = '*** WARNING: unable to remove '+filename
        #print(msg)
        raise OSError(msg)
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise OSError(msg) # re-raise exception if a different error occurred
    except Exception as e:
        raise e
    return
