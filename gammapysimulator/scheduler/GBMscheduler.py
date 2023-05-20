#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from gammapy.data import GTI
from gammapy.datasets import SpectrumDataset, Datasets
from gammapy.irf import EDispKernel, EDispKernelMap
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
    
    def ScheduleDatasets(self):
        """
        Load IRFs compatible with XSPEC and create a SpectrumDataset with only IRFs.
        """
        
        # Set PSF Containment false
        self.psf_containment = False
        
        # Load GTIs
        GTIs = self.DefineGTIs()
        
        # Load Background Spectrum
        if "BAK" in self.conf.IRFfilepath.keys():
            background=self.ReadBackgroundFromBAK(self.conf.IRFfilepath['BAK'])
        else:
            self.log.warning("Key \'BAK\' not found: no background file selected.")
        
        # Load IRFs: Effectiva Area and Energy Dispersion Probability
        if "RSP" in self.conf.IRFfilepath.keys():
            DetectorResponseMatrix = self.ReadDRMfromRSP(self.conf.IRFfilepath['RSP'])
            exposure= self.ComputeExposureFromDRM(DetectorResponseMatrix)
            edisp   = self.ComputeEDispProFromDRM(DetectorResponseMatrix)
        else:
            exposure= self._read_Exposure_from_ARF(self.conf.IRFfilepath['ARF'])
            edisp   = self._read_EDispPro_from_RMF(self.conf.IRFfilepath['RMF'])
        
        # Define a collection of empty Datasets with IRFs and GTIs.
        self.emptydatasets = self.SetEmptyDatasets(GTIs, background, exposure, edisp)
        
        return None
        

    def DefineGTIs(self):
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
    
    def ReadBackgroundFromBAK(self, Filepath):
        """
        Read Background Spectrum from a .bak file.
        
        Parameters
        ----------
        Filepath : pathlib.Path
            Path of the background file.
        
        Return
        ------
        CountsMap : gammapy.maps.RegionNDMap
            Map containing the background counts spectrum.
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
                                               unit=self.conf.energyUnit,
                                               name='energy'
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
        RatesInterpolated = utils.InterpolateFunction(Rates,
                                                      EnergyAxis.center,
                                                      self.conf.AxisEnergyReco.center,
                                                      scale='log'
                                                      )
        
        # Plot the old and new Rates
        self.exporter.PlotStep([(EnergyAxis.center, Rates, "Rates File"),
                              (self.conf.AxisEnergyReco.center, RatesInterpolated, "Rates Interpolated")
                              ],
                             {'xlabel': f"Energy / {self.conf.AxisEnergyReco.unit}",
                              'ylabel': f"Rates / ({RatesInterpolated.unit})",
                              'xscale' : "log",
                              'yscale' : "log",
                              'title' : "Background Rates IRF",
                              'figurename' : "irfs/BackgroundRates_bak"
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
                              'figurename' : "irfs/BackgroundCountsSpectrum"
                              }
                             )     
        
        return CountsMap
    
    
    
    def ReadDRMfromRSP(self, Filepath):
        """
        Read Background Spectrum from a .bak file.
        
        Parameters
        ----------
        Filepath : pathlib.Path
            Path of the background file.
        
        Return
        ------
        MapDRM : gammapy.maps.RegionNDMap
            Map containing the Detector Response Matrix and its geometry.
        """
        # Check file existence
        self.log.info(f"Load DRM from RSP file: {Filepath}")
        if not os.path.isfile(Filepath):
            raise FileNotFoundError(f"File Not Found: {Filepath}")    
        
        # Load requested HDUs
        with fits.open(Filepath) as hdulist:
            Ebounds = Table.read(hdulist['EBOUNDS' ])
            SpecResp= Table.read(hdulist['SPECRESP MATRIX'])

        # Create the Reco Energy Axis of the Matrix
        EnergyAxisReco = MapAxis.from_energy_edges(np.append(Ebounds['E_MIN'],
                                                             Ebounds['E_MAX'][-1]
                                                             ),
                                                   unit=Ebounds['E_MIN'].unit,
                                                   name='energy'
                                                   )
        
        # Create the True Energy Axis of the Matrix
        EnergyAxisTrue = MapAxis.from_energy_edges(np.append(SpecResp['ENERG_LO'],
                                                             SpecResp['ENERG_HI'][-1]
                                                             ),
                                                   unit=SpecResp['ENERG_LO'].unit,
                                                   name='energy_true'
                                                   )
        
        # Define the Geometry
        Geometry = RegionGeom.create(CircleSkyRegion(self.conf.target, self.conf.RegionRadius),
                                     axes = [EnergyAxisReco, EnergyAxisTrue]
                                     )
        
        # Load the 2D Matrix
        DRM = np.zeros([EnergyAxisTrue.nbin, EnergyAxisReco.nbin], dtype = np.float64)
    
        for i, rowSpecResp in enumerate(SpecResp):
            if rowSpecResp["N_GRP"]:
                m_start = 0
                for k in range(rowSpecResp["N_GRP"]):
                
                    if np.isscalar(rowSpecResp["N_CHAN"]):
                        f_chan = rowSpecResp["F_CHAN"]    -1 # Necessary only for GBM (?)
                        n_chan = rowSpecResp["N_CHAN"]
                    else:
                        f_chan = rowSpecResp["F_CHAN"][k] -1 # Necessary only for GBM (?)
                        n_chan = rowSpecResp["N_CHAN"][k]

                    DRM[i, f_chan : f_chan+n_chan] = rowSpecResp["MATRIX"][m_start : m_start+n_chan]
                    m_start += n_chan

        # Add Unit to the DRM Matrix
        DRMUnit = u.Unit(SpecResp['MATRIX'].unit)

        # Combine Geometry and DRM Matrix into RegionNDMap
        MapDRM = RegionNDMap.from_geom(Geometry, data=DRM, unit=DRMUnit, meta=SpecResp.meta)
        
        try:
            MapDRM = MapDRM.to_unit('m2')
        except:
            self.log.warning(f"DRM [{DRMUnit}] conversion to m2 failed")
        
        self.exporter.PlotDRM(MapDRM, stretch=0.3, filename="DetectorResponseMatrix_rsp")
        
        return MapDRM
    
    
    def ComputeExposureFromDRM(self, MapDRM):
        """
        Compute an Exposure Map in the analysis Geometry from the MapDRM.
        
        Parameters
        ----------
        MapDRM : gammapy.maps.RegionNDMap
            Map containing the Detector Response Matrix and its geometry.
        
        Return
        ------
        ExposureMap : gammapy.maps.RegionNDMap
            Map containing the exposure.
        """
        # Compute Effective Area by summing DRM over reco energy axis
        Aeff = np.sum(np.squeeze(MapDRM.data), axis=1)
        Aeff = Aeff * MapDRM.unit
        
        # This Map has a different precision wrt to the requested analysis binning,
        # we must interpolate rates.
        AeffInterpolated = utils.InterpolateFunction(Aeff,
                                                     MapDRM.geom.axes['energy_true'].center,
                                                     self.conf.AxisEnergyTrue.center,
                                                     scale='log'
                                                     )
        
        # Plot the old and new Rates
        self.exporter.PlotStep([(MapDRM.geom.axes['energy_true'].center, Aeff, "AEFF File"),
                                (self.conf.AxisEnergyTrue.center, AeffInterpolated, "AEFF Interpolated")
                                ],
                               {'xlabel': f"Energy True / {self.conf.AxisEnergyTrue.unit}",
                                'ylabel': f"Effective Area / ({AeffInterpolated.unit})",
                                'xscale': "log",
                                'yscale': "log",
                                'title' : "Effective Area IRF",
                                'figurename' : "irfs/EffectiveArea_rsp"
                                }
                               )
        
        # Get the Exposure by Multiplying Effective Area for the simulation integration time.
        Exposure = AeffInterpolated * self.conf.timeReso
        
        # Define the Exposure Geometry
        Geometry = RegionGeom.create(CircleSkyRegion(self.conf.target, self.conf.RegionRadius),
                                     axes = [self.conf.AxisEnergyTrue]
                                     )
        
        # Create the RegionNDMap with Counts and plot
        ExposureMap = RegionNDMap.from_geom(Geometry,
                                            data=Exposure.value,
                                            unit=Exposure.unit,
                                            meta=MapDRM.meta
                                            )
        ExposureMap.meta['livetime'] = self.conf.timeReso
        self.exporter.PlotStep([(self.conf.AxisEnergyTrue.center, Exposure, "Exposure")],
                               {'xlabel': f"Energy True / {self.conf.AxisEnergyTrue.unit}",
                                'ylabel': f"Exposure / {Exposure.unit}",
                                'xscale' : "log",
                                'yscale' : "log",
                                'title' : f"Exposure IRF [{self.conf.timeReso}]",
                                'figurename' : "irfs/Exposure"
                                }
                               )
        
        return ExposureMap
    
    def ComputeEDispProFromDRM(self, MapDRM):
        """
        Compute Energy Dispersion Matrix from a Detector Response Matrix.
        
        Parameters
        ----------
        MapDRM : gammapy.maps.RegionNDMap
            Map containing the Detector Response Matrix and its geometry.
        
        Return
        ------
        EnergyDispersionMap : `gammapy.irf.EDispKernelMap `
            Energy Dispersion Matrix.
        """
        # MapDRM has a different resolution than the one requested for the simulation.
        # Perform 2D interpolation of the Effective Area Matrix.
        DRMInterpolated = utils.InterpolateMap(OldValues=np.squeeze(MapDRM.data)* MapDRM.unit,
                                               OldAxis1=MapDRM.geom.axes['energy_true'].center,
                                               OldAxis2=MapDRM.geom.axes['energy'].center,
                                               NewAxis1=self.conf.AxisEnergyTrue.center,
                                               NewAxis2=self.conf.AxisEnergyReco.center,
                                               scale1='log',
                                               scale2='log'
                                               )
        DRM = DRMInterpolated.value
        DRMUnit = DRMInterpolated.unit
        
        # Define Geometry with two energy axes
        Geometry = RegionGeom.create(CircleSkyRegion(self.conf.target, self.conf.RegionRadius),
                                     axes = [self.conf.AxisEnergyReco, self.conf.AxisEnergyTrue]
                                     )

        # Define full exposure map
        # Compute Effective Area by summing DRM over reco energy axis
        exposure = np.sum(DRM, axis=1) * DRMUnit * self.conf.timeReso
        
        exposure_map = RegionNDMap.from_geom(Geometry.squash('energy'),
                                             data=exposure.value,
                                             unit=exposure.unit,
                                             meta=MapDRM.meta
                                             )
        
        # Normalize Detector Response Matrix to get a probability distribution
        self.log.warning("Normalizing DRM columns to 1: assumption of no photons lost due to energy dispersion outside energy thresholds.")
        ChannelNorm = np.sum(DRM, axis=1)
        for i, norm in enumerate(ChannelNorm):
            if norm > 0.0:
                DRM[i] = DRM[i] / norm
        
        # EDispKernelMap from constructor
        edisp_map = RegionNDMap.from_geom(Geometry, data=DRM, unit='', meta=MapDRM.meta)
        EnergyDispersionMap = EDispKernelMap(edisp_kernel_map=edisp_map, exposure_map=exposure_map)
        
        # Plot the Energy Dispersion
        self.exporter.PlotDRM(EnergyDispersionMap.edisp_map, stretch=0.3, cmap='viridis', filename="EnergyDispersionProbability_rsp")
        
        return EnergyDispersionMap
    
    def _read_Exposure_from_ARF(self, filepath):
        raise NotImplementedError
    
    def _read_EDispPro_from_RMF(self, filepath):
        raise NotImplementedError
    
    def SetEmptyDatasets(self, GTIs, BackgroundMap, ExposureMap, EnergyDispersionMap):
        """
        Define a collection of SpectrumDatasets with different GTIs but same DL4 IRFs.
        
        Parameters
        ----------
        GTIs : list of gammapy.data.GTI
            GTIs that define the time bin of each observation.
        BackgroundMap : gammapy.maps.RegionNDMap
            Map containing the background counts spectrum.
        ExposureMap : gammapy.maps.RegionNDMap
            Map containing the exposure.
        EnergyDisperionMap : gammapy.maps.RegionNDMap
            Map containing the energy dispersion probability.
            
        Returns
        -------
        datasets :  gammapy.datasets.Datasets()
            A collection of SpectrumDataset with no counts.
        """
        datasets = Datasets()
        EmptyCounts = RegionNDMap.from_geom(BackgroundMap.geom,
                                            data=0,
                                            unit='',
                                            )
        MaskSafe = RegionNDMap.from_geom(BackgroundMap.geom,
                                         data=True,
                                         )
        for i, gti in enumerate(GTIs):
            spectrumdataset = SpectrumDataset(counts=EmptyCounts,
                                              background=BackgroundMap,
                                              exposure=ExposureMap,
                                              edisp=EnergyDispersionMap,
                                              gti=gti,
                                              name=f"dataset-{i}",
                                              mask_safe=MaskSafe
                                              )
            datasets.append(spectrumdataset)
        
        return datasets
