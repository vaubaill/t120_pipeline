import logging
from astropy.time import Time

def get_module_logger(mod_name,logfile='',addsldr=True):
    """Creates a logger object.
    
    Parameters
    ----------
    mod_name : string
        module name.
    logfile : string, optional
        Log file name. Default is provided by the get_log_file_name() function.
    addsldr : Booleanm optional
        If True, a streamHandler is added to the logger object.
        Default is True.
    
    Returns
    -------
    log : logging.logger object, with file handler and optionaly stream handler.
    
    """
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
    """Creates a logger file name from current date.
    
    Parameters
    ----------
    pathlog : string, optional
        Log file directory name.
    
    Returns
    -------
    log_file : string
        Log file full name.
    
    """
    log_file = pathlog + 't120-'+Time.now().isot.replace('-','').replace(':','').split('.')[0]+'.log'
    return log_file
