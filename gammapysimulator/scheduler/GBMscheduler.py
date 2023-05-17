#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

from astropy import units as u
from astropy.io import fits
from astropy.table import Table, QTable
from gammapy.data import GTI
from gammapy.maps import RegionNDMap, MapAxis, RegionGeom
from regions import CircleSkyRegion

import numpy as np
import os

from gammapysimulator.configure.configure import SimulationConfigurator
from gammapysimulator.scheduler.scheduler import Scheduler
from gammapysimulator.tools.export import ExportSimulations
from gammapysimulator.tools import utils


class GBMScheduler(Scheduler):
    """
    Class that organises Observations and Datasets with Fermi GBM IRFs.
    """
    
    def __init__(self, configurator : SimulationConfigurator, exporter : ExportSimulations) -> None:
        """
        Instantiate Scheduler by reading IRFs.
        
        Parameters
        ----------
        configurator : gammapysimulator.configure.configure.SimulationConfigurator
            Configurator object.
        exporter : gammapysimulator.tools.export.ExportSimulations
            Export object.
        """
        
        self.conf= configurator
        self.log = configurator.log
        self.exporter = exporter
    
    def LoadIRFs(self):
        """
        Load IRFs compatible with XSPEC and create a SpectrumDataset with only IRFs.
        """
        
        # Load GTIs
        GTIs = self._define_GTIs()
        
        # Load Background Spectrum
        if "BAK" in self.conf.IRFfilepath.keys():
            background=self._read_Background_from_BAK(self.conf.IRFfilepath['BAK'])
        else:
            self.log.warning("Key \'BAK\' not found: no background file selected.")
        
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
        self.emptydatasets = self._Set_Empty_Datasets()
        
        return None
        

    def _define_GTIs(self):
        """
        Compute the Schedule and define the GTI objects for the simulations.
        Each one represents a different time bin of the simulation.
        
        Return
        ------
        GTIs : list of gammapy.data.GTI
            List of GTI objects.
        """
        ObservationsStart, ObservationsStop = self.DefineSchedule()
        GTIs = [GTI.create(start, stop, reference_time=self.conf.timeRef) for start, stop in zip(ObservationsStart, ObservationsStop)]
        
        return GTIs
    
    def _read_Background_from_BAK(self, Filepath):
        """
        Read Background Spectrum from a .bak file.
        
        Parameters
        ----------
        Filepath : pathlib.Path
            Path of the background file.
        
        Return
        ------
        background : gammapy.maps.RegionNDMap
            Map containing the background spectrum.
        """
        # Check file existence
        self.log.info(f"Load BAK from file: {Filepath}")
        if not os.path.isfile(Filepath):
            raise FileNotFoundError(f"File Not Found: {Filepath}")    
        
        # Load requested HDUs
        with fits.open(Filepath) as hdulist:
            Ebounds = Table.read(hdulist['EBOUNDS' ])
            Spectrum= Table.read(hdulist['SPECTRUM'])

        # Create the energy axis of the Spectrum and the Geometry
        EnergyAxis = MapAxis.from_energy_edges(np.append(Ebounds['E_MIN'],
                                                         Ebounds['E_MAX'][-1]
                                                         ),
                                               unit=self.conf.energyUnit
                                               )
        Geometry = RegionGeom.create(CircleSkyRegion(self.conf.target, self.conf.RegionRadius),
                                     axes = [EnergyAxis]
                                     )
        
        # Get the background rates, and scale them according to the simulation integration time
        if 'RATE' in Spectrum.colnames:
            RateUnit = u.Unit(f"{Spectrum.meta['TIMEUNIT']}-1")
            Rates = Spectrum['RATE'] * RateUnit
        elif 'COUNTS' in Spectrum.colnames:
            IntegrationTime = 1.0*u.s
            Rates = Spectrum['COUNTS'].data / IntegrationTime
            raise NotImplementedError("READ INTEGRATION TIME OF THE COUNTS")
        else:
            raise KeyError()
        
        # This Map has a different precision wrt to the requested analysis binning,
        # we must interpolate rates.
        RatesInterpolated = utils.InterpolateMap(Rates, EnergyAxis.center, self.conf.AxisEnergyReco.center, scale='log')
        
        # Plot the old and new Rates
        self.exporter.PlotStep([(EnergyAxis.center, Rates, "Rates File"),
                              (self.conf.AxisEnergyReco.center, RatesInterpolated, "Rates Interpolated")
                              ],
                             {'xlabel': f"Energy / {self.conf.AxisEnergyReco.unit}",
                              'ylabel': f"Rates / ({RatesInterpolated.unit})",
                              'xscale' : "log",
                              'yscale' : "log",
                              'title' : "Background Rates IRF",
                              'figurename' : "irfs/ratesbak"
                              }
                             )
        
        # Get the Counts by Multiplying Rates for the simulation integration time.
        Counts = (RatesInterpolated * self.conf.timeReso).to('') 
        
        # Create the RegionNDMap with Counts and plot
        CountsMap = RegionNDMap.from_geom(self.conf.geometry, data=Counts.value, unit=Counts.unit, meta=Spectrum.meta)
        self.exporter.PlotStep([(self.conf.AxisEnergyReco.center, Counts, "Background")],
                             {'xlabel': f"Energy / {self.conf.AxisEnergyReco.unit}",
                              'ylabel': f"Counts",
                              'xscale' : "log",
                              'yscale' : "log",
                              'title' : f"Background Counts IRF [{self.conf.timeReso}]",
                              'figurename' : "irfs/countsbak"
                              }
                             )     
        
        return CountsMap
    
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
