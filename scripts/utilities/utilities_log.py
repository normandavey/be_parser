import logging,os,inspect

file_path = os.path.dirname(inspect.stack()[0][1])

###################################################################
### LOGS
###################################################################

def setupLogger(options):
    if 'debug' not in options: options['debug'] = False

    if options['debug']:
        logging_level = logging.DEBUG
    elif options['verbose']:
        logging_level = logging.INFO
    else:
        logging_level = logging.ERROR

    if 'logfile' not in options: options['logfile'] = False
    elif options['logfile'] == "": options['logfile'] = False
    elif isinstance(options['logfile'],(str)): pass
    elif options['logfile']: options['logfile'] = os.path.join(file_path,'data_manager.log')
    
    # create console handler with a higher log level
    formatter = logging.Formatter('%(asctime)s %(levelname)-8s [%(funcName)s] %(message)s')
    ch = logging.StreamHandler()
    ch.setLevel(logging_level)
    ch.setFormatter(formatter)

    if options['logfile'] != False:
        # create file handler which logs even debug messages
        formatter = logging.Formatter('%(asctime)s %(levelname)-8s [%(funcName)s] %(message)s')
        fh = logging.FileHandler(options['logfile'])
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)

        logging.basicConfig(
            format='%(asctime)s %(levelname)-8s [%(funcName)s] %(message)s',
            datefmt='%H:%M:%S',
            level=logging.DEBUG,
            handlers=[ch,fh]
            )
        
        logger = logging.getLogger(__name__)
        logger.info("Writing log: " + options['logfile'])
    
    else:
        logging.basicConfig(
            format='%(asctime)s %(levelname)-8s [%(funcName)s] %(message)s',
            datefmt='%H:%M:%S',
            level=logging_level
            )

        logger = logging.getLogger(__name__)

    return __name__ 