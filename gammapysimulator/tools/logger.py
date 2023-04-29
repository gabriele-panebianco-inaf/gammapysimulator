#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

import logging
from pathlib import Path
from time import strftime

class Singleton(type):
    '''Make sure there is a single instance of the logger class at any time.'''
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(
                Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]

class SimulatorLogger(metaclass=Singleton):
    """Define a class for logging"""
    
    def __init__(self,
                 LoggerName="simulator",
                 LogLevel='DEBUG',
                 LogToFile=False,
                 LogOutDir='log',
                 LogFormat='%(asctime)s - %(name)s - %(levelname)s - %(message)s'                 
                 )->None:      
        """
        Initialise logger defining logging level, format and output.
        
        Parameters
        ----------
        LoggerName : str
            Logger Name.
        LogLevel : str
            Logging level.
        LogToFile :  bool
            Decide if to log on a file.
        LogOutDir : str
            path to destination folder.
        """        
        
        self.loggerName = LoggerName
        self.formatter = logging.Formatter(LogFormat)
        self.logLevel = LogLevel
        
        self.logToFile = LogToFile
        if self.logToFile:
            self.logOutDir = LogOutDir
            Path(LogOutDir).mkdir(parents=True, exist_ok=True)

        self.streamHandler = logging.StreamHandler()
        self.streamHandler.setLevel(self.logLevel)
        self.streamHandler.setFormatter(self.formatter)

    def getLogger(self):
        """Getter for the Logger.

        Return
        ------
        logger : SimulatorLogger
            Logger.
        """        
        
        logger = logging.getLogger(self.loggerName)
        
        logger.setLevel(self.logLevel)
        if not logger.handlers:
            logger.addHandler(self.streamHandler)
        logger.propagate = False

        if self.logToFile:
            self.fh = logging.FileHandler(str(Path(self.logOutDir).joinpath(f"{self.loggerName}_{strftime('%Y%m%d-%H%M%S')}.log")))
            self.fh.setLevel(self.logLevel)
            self.fh.setFormatter(self.formatter)
            if not logger.handlers:
                logger.addHandler(self.fh)
        
        return logger

    