#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

from astropy import units as u
from astropy.table import Table
from astropy.time import Time
from pathlib import Path
from time import strftime

import numpy as np

from gammapysimulator.tools import logger

class FormatLightCurve():
    """
    Class that reads a light curve from different formats and write them as
    a fits file that can be read by gammapy as a gammapy.modeling.models.LightCurveTemplateTemporalModel.
    """
    def __init__(self, lightcurvefile, inputformat, referenceMJD) -> None:
        """
        Read the light curve file in an astropy.table.Table
        
        Parameters
        ----------
        lightcurvefile : str
            Path to light curve file.
        inputformat : str
            Format of the input light curve.
        referenceMJD : str
            MJD reference of the light curve.
        """
        # Configure logger
        self.log = logger.SimulatorLogger().getLogger()
        
        # Set inputfile
        self.inputfile = Path(lightcurvefile)
        if not self.inputfile.is_file():
            raise FileNotFoundError
        
        # Set MJD referencetime
        self.log.info(f"Using reference time MJD={referenceMJD}")
        self.TimeRef = Time(referenceMJD, format='mjd')
        
        # Read Table
        if inputformat=="megalib":
            self.table = self.ReadFromMegaLib()
        else:
            raise NotImplementedError("Only megalib format currently supported")
        
        return None
    
    def write(self, outputfilename, overwrite=True):
        """Write table to an output destination.
        
        Parameters
        ----------
        outputfilename : str
            Path of the output light curve file.
        overwrite : bool
            Overwrite parameter.
        """
        self.log.info(f"Write {outputfilename}")
        self.table.write(outputfilename, format='fits', overwrite=overwrite)
        return None
    
    def ReadFromMegaLib(self):
        """
        Read a light curve file written in a format compatible with MegaLib
            
        Return
        ------
        table : astropy.table.Table
            Table with the formatted light curve.
        """
        self.log.info("Read a light curve written in MegaLib format.")
        
        self.log.info(f"Read light curve from {self.inputfile}")
        
        input_table = Table.read(self.inputfile,
                                 guess=False,
                                 fast_reader=False,
                                 format='ascii',
                                 comment='#',
                                 delimiter='\s',
                                 header_start=None,
                                 data_start=1,
                                 data_end=-1,
                                 names = ["DP","time","counts"]
                                 )
        
        # Define output table
        table = Table()
        
        # Set TIME column
        table['TIME']=input_table['time']
        
        # Set NORM column, that must contain values in [0,1]
        MaxCounts = np.max(input_table['counts'])
        table['NORM']=input_table['counts']/MaxCounts
        
        # Set Metadata
        table.meta['EXTNAME']="Time profile"
        table.meta['MJDREFI']= np.floor(self.TimeRef.value)
        table.meta['MJDREFF']= self.TimeRef.value - table.meta['MJDREFI']
        table.meta['TIMEUNIT']="s"
        table.meta['TIMESYS']="TT"
        table.meta['TIMEREF']="LOCAL"
        table.meta['MAXNORM']=MaxCounts
        
        return table