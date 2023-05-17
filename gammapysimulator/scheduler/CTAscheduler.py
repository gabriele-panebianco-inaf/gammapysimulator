#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

from gammapy.data import Observation, Observations
from gammapy.datasets import SpectrumDataset, MapDataset, Datasets
from gammapy.irf import load_cta_irfs
from gammapy.makers import SpectrumDatasetMaker, MapDatasetMaker, SafeMaskMaker
from time import time
from tqdm import tqdm

from gammapysimulator.configure.configure import SimulationConfigurator
from gammapysimulator.scheduler.scheduler import Scheduler
from gammapysimulator.tools.export import ExportSimulations

import numpy as np

class CTAScheduler(Scheduler):
    """
    Class that organises Observations and Datasets with CTA IRFs.
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
        
        # Load IRFs
        self.log.info(f"Load CTA IRFs from file: {self.conf.IRFfilepath}")      
        self.irfs = load_cta_irfs(self.conf.IRFfilepath)
        self.exporter.PlotCTAIRFs(self.irfs)

    def ScheduleDatasets(self):
        """
        Load IRFs compatible with XSPEC and create a SpectrumDataset with only IRFs.
        """
        # Read IRFs and Set the observations
        observations = self.SetObservations()
        
        # Set PSF Containment
        self.psf_containment = "psf" in self.irfs
        self.log.info(f"PSF Containment: {self.psf_containment}")
        
        # Reduce IRFs to get the Datasets
        self.emptydatasets = self.MakeReducedIRFsDatasets(observations)
        
        return None


    def SetObservations(self):
        """
        Define the Observations object.
        
        Return
        ------
        observations : gammapy.data.Observations
            Collection of Observation with IRFs and time bin information.
        """
        
        # Create Schedule
        ObservationsStart, ObservationsStop = self.DefineSchedule()
        
        # Define observations
        observations = Observations()
        
        # Observations are time bins, they differ by ID, start and stop
        for i, (tstart, tstop) in enumerate(zip(ObservationsStart, ObservationsStop)):
            obs = Observation.create(self.conf.pointing,
                                     obs_id=i,
                                     tstart=tstart,
                                     tstop=tstop,
                                     irfs=self.irfs,
                                     reference_time=self.conf.timeRef
                                     )
            observations.append(obs)
            
        return observations

    def MakeReducedIRFsDatasets(self, observations):
        """
        Bundle IRFs and GTIs into a collection of Datasets.
        
        Parameters
        ----------
        observations : gammapy.data.Observations
            Collection of observations with DL3 IRFs and GTIs.
        
        Return
        ------
        datasets : gammapy.datasets.Datasets
            Collection of Datasets with Reduced IRFs and GTIs.
        """
        if self.conf.analysis=="3D":
            empty = MapDataset.create(geom = self.conf.geometry,
                                      energy_axis_true = self.conf.AxisEnergyTrue,
                                      name = "empty",
                                      )
            maker = MapDatasetMaker(selection = ['exposure', 'background', 'psf', 'edisp'])
        
            safe_mask_maker = SafeMaskMaker(methods = ["aeff-default","offset-max"],
                                            offset_max = 2*self.conf.FoVRadius
                                            )
        elif self.conf.analysis=="1D":
            empty = SpectrumDataset.create(geom = self.conf.geometry,
                                           energy_axis_true = self.conf.AxisEnergyTrue,
                                           name = "empty"
                                           )
            maker = SpectrumDatasetMaker(selection = ["exposure", "background", "edisp"],
                                         containment_correction = self.psf_containment
                                         )
            safe_mask_maker = SafeMaskMaker(methods = ["aeff-default"])
        
        # Define Datasets collection
        datasets = Datasets()
        
        # Perform reduction and measure time
        self.log.info(f"Perform DL3->DL4 IRFs reduction for a {self.conf.analysis} simulation...")
        now = time()
        
        for idx, obs in enumerate(tqdm(observations, desc="Reduction Loop")):

            # Set Name and run DL3->DL4 Reduction
            dataset = maker.run(empty.copy(name = f"dataset-{idx}"), obs)
            dataset = safe_mask_maker.run(dataset, obs)
            
            # Add current dataset to the Datasets() object
            datasets.append(dataset)
            
        self.log.info(f"Reduction performed in {float(time()-now):.3f} s.")
        
        # Return the simulated Datasets
        return datasets
        