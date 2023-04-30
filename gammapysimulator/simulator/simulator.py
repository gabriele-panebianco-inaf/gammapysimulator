#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

from gammapy.datasets import MapDataset, SpectrumDataset, Datasets
from gammapy.makers import MapDatasetMaker, SpectrumDatasetMaker, SafeMaskMaker
from gammapy.modeling.models import PowerLawSpectralModel, PointSpatialModel, SkyModel, ConstantTemporalModel
from tqdm import tqdm
from time import time

from gammapysimulator.configure.configure import SimulationConfigurator
from gammapysimulator.scheduler import CTAscheduler


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
        configurator = SimulationConfigurator()
        configurator.read(ConfigurationFileName)
        self.conf = configurator
        
        # Set Log
        self.log = self.conf.log
        
        # Set Observations
        self.observations = self.SetObservations()
        
        # Set the Models
        self.models = self.SetModels()
        

    def SetObservations(self):
        """
        Set the Observations object according to the Instrument.
        """
        
        if self.conf.instrument=="CTA":
            scheduler = CTAscheduler.CTAScheduler(self.conf)
            scheduler.SetObservations()
            observations = scheduler.observations
        else:
            raise NotImplementedError
        
        # Set PSF Containment
        self.psf_containment = "psf" in scheduler.irfs
        self.log.info(f"PSF Containment: {self.psf_containment}")
        
        return observations
    
    def SetModels(self):
        """
        Set the Simulation Models.
        """
        spectral = PowerLawSpectralModel(index=2.4,
                                         amplitude="5.7e-11 cm-2 s-1 TeV-1",
                                         reference="1 TeV"
                                         )
        
        temporal = ConstantTemporalModel(const=1)
        
        spatial = PointSpatialModel(lon_0="83.00 deg", lat_0="22.50 deg",frame="fk5")
        
        models = SkyModel(spectral_model = spectral,
                          spatial_model = spatial,
                          temporal_model=temporal,
                          name = "Current")

        return models
        
    def SetDatasets(self):
        """
        Set the Datasets object according to 3D or 1D Analysis.
        """
        
        if self.conf.analysis=="3D":
            empty = MapDataset.create(geom = self.conf.geometry,
                                      energy_axis_true = self.conf.axis_energy_true,
                                      name = "empty",
                                      )
            maker = MapDatasetMaker(selection = ['exposure', 'background', 'psf', 'edisp'])
        
            safe_mask_maker = SafeMaskMaker(methods = ["aeff-default","offset-max"],
                                            offset_max = 2*self.conf.FoVRadius
                                            )
        elif self.conf.analysis=="1D":
            empty = SpectrumDataset.create(geom = self.conf.geometry,
                                           energy_axis_true = self.conf.axis_energy_true,
                                           name = "empty"
                                           )
            maker = SpectrumDatasetMaker(selection = ["exposure", "background", "edisp"],
                                         containment_correction = self.psf_containment
                                         )
            safe_mask_maker = SafeMaskMaker(methods = ["aeff-default"])
        
        # Perform reduction and measure time
        self.log.info(f"Perform Dataset Reduction...")
        now = time()
        
        datasets = Datasets()
            
        for idx, obs in enumerate(tqdm(self.observations, desc="Reduction Loop")):

            # Set Name and run DL3->DL4 Reduction
            dataset = maker.run(empty.copy(name = f"dataset-{idx}"), obs)
            dataset = safe_mask_maker.run(dataset, obs)

            # Set the source model and compute the predicted excess counts.
            dataset.models = self.models

            datasets.append(dataset)
            
        self.log.info(f"Reduction performed in {float(time()-now):.3f} s.")
        
        # Set Datasets
        self.datasets = datasets
