import logging
from astropy.time import Time

def get_module_logger(mod_name,logfile='',addsldr=True):
    # create logger
    log = logging.getLogger(mod_name)
    # log format
    log_fmt = '%(levelname)s %(filename)s %(lineno)d (%(funcName)s) : %(message)s '
    # set log level
    log.setLevel(logging.DEBUG)
    log.propagate = True
    # create StreamHandler
    sdlr = logging.StreamHandler()
    sdlr.setLevel(logging.DEBUG)
    sdlr.setFormatter(logging.Formatter(fmt=log_fmt))
    if addsldr:
        log.addHandler(sdlr)
    # create FileHandler
    if logfile:
        hdlr = logging.FileHandler(logfile)
        hdlr.setLevel(logging.DEBUG)
        hdlr.setFormatter(logging.Formatter(fmt=log_fmt))
        log.addHandler(hdlr)
    # log the log file name
    log.info('Log file: '+logfile)
    
    return log


def get_log_file_name(pathlog='./'):
    log_file = pathlog + 't120-'+Time.now().isot.replace('-','').replace(':','').split('.')[0]+'.log'
    return log_file
