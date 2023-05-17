#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

from gammapysimulator.scheduler.scheduler import Scheduler

import numpy as np

class GBMScheduler(Scheduler):
    """
    Class that organises Observations and Datasets with Fermi GBM IRFs.
    """
    
    def __init__(self, configurator) -> None:
        """
        Instantiate Scheduler by reading IRFs.
        
        Parameters
        ----------
        configurator : gammapysimulator.configure.configure.SimulationConfigurator
            Configurator object.
        """
        
        self.conf= configurator
        self.log = configurator.log
    
    def LoadIRFs(self):
        """
        Load IRFs compatible with XSPEC and create a SpectrumDataset with only IRFs.
        """
        
        # Load GTIs
        GTIs = self._define_GTIs()
        
        # Load Background Spectrum
        if "BAK" in self.conf.IRFfilepath.keys():
            background=self._read_Background_from_BAK(self.conf.IRFfilepath['BAK'])
        
        # Load IRFs: Effectiva Area and Energy Dispersion Probability
        if "RSP" in self.conf.IRFfilepath.keys():
            exposure= self._read_Exposure_from_RSP(self.conf.IRFfilepath['RSP'])
            edisp   = self._read_EDispPro_from_RSP(self.conf.IRFfilepath['RSP'])
        else:
            exposure= self._read_Exposure_from_ARF(self.conf.IRFfilepath['ARF'])
            edisp   = self._read_EDispPro_from_RMF(self.conf.IRFfilepath['RMF'])
        
        # Define an empty SpectrumDataset with IRFs
        empty = self._define_Empty_Dataset()
        
        # Define a collection of empty Datasets with IRFs and GTIs.
        self.empties = self._Set_Empty_Datasets()
        
        return None
        

    def _define_GTIs(self):
        ObservationsStart, ObservationsStop = self.DefineSchedule()
        raise NotImplementedError
    
    def _read_Background_from_BAK(self, filepath):
        raise NotImplementedError
    
    def _read_Exposure_from_RSP(self, filepath):
        raise NotImplementedError
    
    def _read_EDispPro_from_RSP(self, filepath):
        raise NotImplementedError
    
    def _read_Exposure_from_ARF(self, filepath):
        raise NotImplementedError
    
    def _read_EDispPro_from_RMF(self, filepath):
        raise NotImplementedError
    
    def _define_Empty_Dataset(self):
        raise NotImplementedError
    
    def _Set_Empty_Datasets(self):
        raise NotImplementedError
