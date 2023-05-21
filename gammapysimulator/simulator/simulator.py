#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

from gammapy.datasets import MapDataset, SpectrumDataset, Datasets, SpectrumDatasetOnOff, MapDatasetOnOff
from gammapy.makers import MapDatasetMaker, SpectrumDatasetMaker, SafeMaskMaker
from gammapy.modeling.models import Models, SkyModel
from numpy.random import RandomState
from tqdm import tqdm
from time import time

from gammapysimulator.configure import configure
from gammapysimulator.scheduler import CTAscheduler, GBMscheduler
from gammapysimulator.tools import export


class Simulator:
    """
    Class that execute simulation of gamma-ray sources with gammapy.
    """
    
    def __init__(self, ConfigurationFileName) -> None:
        """
        Instantiate Simulator by reading a configuration file.
        
        Parameters
        ----------
        ConfigurationFileName : str
            Name of the input YAML file.
        """
        
        # Read Configuration file and Set Configurator
        configurator = configure.SimulationConfigurator()
        configurator.read(ConfigurationFileName)
        self.conf = configurator
        
        # Set Log
        self.log = self.conf.log
        
        # Set Export object
        self.exporter = export.ExportSimulations(self.conf)
        
        # Set Empty Datasets with Reduced IRFs and Time binning information
        self.emptydatasets = self.SetDL4IRFs()
        
        # Set the Models
        self.models = self.SetModels()
        
        return None

    def SetDL4IRFs(self):
        """
        Set the Observations object according to the Instrument.
        
        Return
        ------
        emptydatasets : gammapy.datasets.Datasets
            Collection of Datasets containing temporal and IRFs information.
        """
        
        # Define the Scheduler for the requested instrument
        if self.conf.instrument=="CTA":
            scheduler = CTAscheduler.CTAScheduler(self.conf, self.exporter)
        elif self.conf.instrument=="Fermi-GBM":
            scheduler = GBMscheduler.GBMScheduler(self.conf, self.exporter)
        else:
            raise NotImplementedError
        
        # Define Datasets: read IRFs, GTIs and merge them into
        # a collection of Datasets with only reduced IRFs,
        # ready to apply models and perform simulations.
        scheduler.ScheduleDatasets()
        
        # Get the Datasets with reduced IRFs
        emptydatasets = scheduler.getEmptyDatasets()
        
        # Set PSF Containment
        self.psf_containment = scheduler.psf_containment
        self.log.info(f"PSF Containment: {self.psf_containment}")
        
        return emptydatasets
    
    def SetModels(self):
        """
        Set the Simulation Models.
        
        Return
        ------
        models :  gammapy.modeling.Models
            Models of the Simulation.
        """
        # Read Models from a YAML file
        self.log.info(f"Load Models from: {self.conf.modelfilepath}")
        models = Models.read(self.conf.modelfilepath)
        
        # Print
        self.log.info(models)
        
        # plot
        self.exporter.PlotTemporalModel(models[0].temporal_model)
        
        return models
    
    def RunSimulation(self):
        """
        Run the simulation with CTA IRFs.
        
        Return
        ------
        OnOffDatasets : gammapy.datasets.Datasets()
            Collection of SpectrumDatasetOnOff (1D) or MapDatasetOnOff (3D) with simulated data.
        """
        # Measure time
        self.log.info(f"Perform a {self.conf.analysis} Simulation...")
        now = time()
        
        # Set Random State
        random_state = RandomState(self.conf.seed)
        
        # Define Datasets collection
        EmptyDatasets = self.emptydatasets
        NewDatasets = Datasets()
        
        for idx, dataset in enumerate(tqdm(EmptyDatasets, desc="Simulation Loop")):

            # Set the source model and compute the predicted excess counts.
            dataset.models = self.models

            # Fake ON counts
            dataset.fake(random_state=random_state)
            
            # Fake OFF counts realization
            #if self.conf.analysis=="3D":
            #    dataset_onoff = MapDatasetOnOff.from_map_dataset(dataset=dataset,
            #                                                     acceptance=1,
            #                                                     acceptance_off=5,
            #                                                     name = f"onoff-{idx}"
            #                                                     )
            #elif self.conf.analysis=="1D":
            #    dataset_onoff = SpectrumDatasetOnOff.from_spectrum_dataset(dataset=dataset,
            #                                                               acceptance=1,
            #                                                               acceptance_off=5,
            #                                                               name = f"onoff-{idx}"
            #                                                               )
            #dataset_onoff.fake(npred_background=dataset.npred_background(),
            #                   random_state=random_state
            #                   )
            
            # Add current dataset to the Datasets() object
            NewDatasets.append(dataset)
            
        self.log.info(f"Reduction performed in {float(time()-now):.3f} s.")
        
        # Return the simulated Datasets
        return NewDatasets
