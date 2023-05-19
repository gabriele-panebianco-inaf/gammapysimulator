#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from gammapy.maps import RegionGeom, WcsGeom, MapAxis
from pathlib import Path
from regions import CircleSkyRegion

import logging
import yaml
import shutil

from gammapysimulator.tools import logger, utils

AllowedProducts = ["DL3","DL4"]
AllowedAnalysis = ["1D" ,"3D" ]
AllowedInstrument = ["CTA", "Fermi-GBM"]

class SimulationConfigurator:
    """
    Class that reads a YML file and defines the parameters needed for the simulation.
    """
    
    def __init__(self)->None:
        
        self.log = logger.SimulatorLogger().getLogger()

    
    def read(self, ConfigurationFileName):
        """
        Read the YAML Configuration file into a dictionary.
        Create the output directory and a log file.
    
        Parameters
        ----------
        ConfigurationFileName : str
            Name of the input YAML file.
        """

        # Load YAML as a dict
        self.log.info(f"Load Configuration YAML: {ConfigurationFileName}")
        with open(ConfigurationFileName) as f:
            configuration = yaml.load(f, Loader = yaml.SafeLoader)
        
        # Set Output Directory and log file
        self.SetOutput(configuration)
        shutil.copyfile(ConfigurationFileName, self.OutputDirectory.joinpath("configuration.yaml"))
            
        # Set General Paramenters of Simulation
        self.SetSimulation(configuration)
        
        #TODO: Develop and Test more cases
        if self.product=="DL3":
            raise NotImplementedError
        
        # Set Model filepath
        self.modelfilepath = Path(utils.get_absolute_path(configuration['Model']['FilePath'])).absolute()
        if not self.modelfilepath.is_file():
            raise FileNotFoundError(f"Model file not found: {self.modelfilepath}")
        
        # Set Instrument and IRFs
        self.SetInstrument(configuration)
        
        # Set Time, Energy and Space Geometry
        self.SetGeometry(configuration)
        
        return None
    
    def SetOutput(self, configuration : dict):
        """
        Set the Output Directory and log to file.
        
        Parameters
        ----------
        configuration : dict
            Configuration dictionary.
        """
        
        # Set Output Directory Name and Tag
        OutputDirectoryPath = Path(utils.get_absolute_path(configuration['Simulation']['OutputDirectory'])).absolute()
        
        self.OutputID = str(configuration['Simulation']['OutputID'])
        self.OutputDirectory = OutputDirectoryPath.joinpath(self.OutputID)
        
        # Output Directory Cleanup
        utils.clean_directory(self.OutputDirectory, self.log)
        
        # Add Log File
        OutputLogFile = str(self.OutputDirectory.joinpath(f"{self.log.name}.log"))
        fh = logging.FileHandler(OutputLogFile)
        fh.setLevel(self.log.level)
        fh.setFormatter(self.log.handlers[0].formatter)
        self.log.addHandler(fh)
        self.log.info(f"Log to {OutputLogFile}")
        
        return None
    
    def SetSimulation(self, configuration : dict):
        """
        Set general parameters for simulation.
        
        Parameters
        ----------
        configuration : dict
            Configuration dictionary.
        """
        
        # Set Seed
        self.seed = int(configuration['Simulation']['RandomSeed'])
        
        # Set Number of Simulations
        self.simN = int(configuration['Simulation']['SimuNumber'])
        
        # Set Product Type
        self.product = str(configuration['Simulation']['Product'])
        if self.product not in AllowedProducts:
            raise ValueError(f"Only {AllowedProducts} allowed, not {self.product}")
        
        # Set Analysis Type
        self.analysis = str(configuration['Simulation']['Analysis'])
        if self.analysis not in AllowedAnalysis:
            raise ValueError(f"Only {AllowedAnalysis} allowed, not {self.analysis}")
        
        return None
    
    def SetInstrument(self, configuration : dict):
        """
        Set Instruemnt and IRFs parameters.
        
        Parameters
        ----------
        configuration : dict
            Configuration dictionary.
        """
        
        # Set Instrument Name
        self.instrument = str(configuration['IRF']['Instrument'])
        if self.instrument not in AllowedInstrument:
            raise ValueError(f"Only {AllowedInstrument} allowed, not {self.instrument}")
        
        # Set Detector Name (Optional)
        try:
            self.detector = str(configuration['IRF']['Detector'])
        except KeyError:
            self.detector = self.instrument
        
        # Check IRF Type and analysis are compatible
        self.irf_pointlike = bool(configuration['IRF']['Pointlike'])
        
        if self.irf_pointlike and (self.analysis=="3D"):
            raise ValueError("Cannot perform 3D simulations with Point-like IRFs.")
        
        # Set IRFs File Path
        try:
            # Single IRFs File: FilePath is a string
            self.IRFfilepath = Path(utils.get_absolute_path(configuration['IRF']['FilePath'])).absolute()
            # Check for file existence
            if not self.IRFfilepath.is_file():
                raise FileNotFoundError(f"IRF file not found: {self.IRFfilepath}")
            
        except TypeError:
            # Multiple IRFs Files (e.g. RMF, ARF, RSP, BAK): FilePath is a dict.
            self.IRFfilepath = {}
            for key in configuration['IRF']['FilePath'].keys():
                self.IRFfilepath[key] = Path(utils.get_absolute_path(configuration['IRF']['FilePath'][key])).absolute()
                # Check for file existence
                if not self.IRFfilepath[key].is_file():
                    raise FileNotFoundError(f"IRF file not found: {self.IRFfilepath}")            
        
        return None
    
    def SetGeometry(self, configuration : dict):
        """
        Set Analysis Time, Space and Energy Geometry.
        
        Parameters
        ----------
        configuration : dict
            Configuration dictionary.
        """
        # Time Parameters
        self.timeUnit = u.Unit(str(configuration['Geometry']['Time']['Unit']))
        self.timeStart= float(configuration['Geometry']['Time']['Start'])* self.timeUnit
        self.timeStop = float(configuration['Geometry']['Time']['Stop']) * self.timeUnit
        self.timeReso = float(configuration['Geometry']['Time']['Resolution'])* self.timeUnit
        self.timeRef  = Time(str(configuration['Geometry']['Time']['Reference']))
        
        # Energy Axes
        self.energyUnit = u.Unit(str(configuration['Geometry']['Energy']['Unit']))
        
        self.AxisEnergyReco = MapAxis.from_energy_bounds(
            configuration['Geometry']['Energy']['RangeReco'][0] * self.energyUnit,
            configuration['Geometry']['Energy']['RangeReco'][1] * self.energyUnit,
            configuration['Geometry']['Energy']['RecoBinPerDecade'],
            per_decade=True,
            name='energy'
            )
        
        self.AxisEnergyTrue = MapAxis.from_energy_bounds(
            configuration['Geometry']['Energy']['RangeTrue'][0] * self.energyUnit,
            configuration['Geometry']['Energy']['RangeTrue'][1] * self.energyUnit,
            configuration['Geometry']['Energy']['TrueBinPerDecade'],
            per_decade=True,
            name='energy_true'
            )
        
        # Space Grid Parameters
        self.frame = str(configuration['Geometry']['Space']['Frame'])
        self.frameUnit = u.Unit(str(configuration['Geometry']['Space']['Unit']))
        
        # Pointing coordinates
        self.pointing = SkyCoord(float(configuration['Pointing']['PointingLon']),
                                 float(configuration['Pointing']['PointingLat']),
                                 unit = self.frameUnit,
                                 frame= self.frame
                                 )
        
        # Target coordinates
        self.target = SkyCoord(float(configuration['Target']['TargetLon']),
                               float(configuration['Target']['TargetLat']),
                               unit = self.frameUnit,
                               frame= self.frame
                               )
        
        # Geometry of the Analysis
        if self.analysis=="3D":
            
            # Center 3D Geometry on Pointing Coordinates
            # Define Spatial Grid by its resolution and FoV width
            
            self.FoVRadius = float(configuration['Geometry']['Space']['FoVRadius']) * self.frameUnit
            self.resolution= float(configuration['Geometry']['Space']['Resolution'])* self.frameUnit

            geometry = WcsGeom.create(skydir= self.pointing,
                                      binsz = self.resolution.to('deg').value,
                                      width = (self.FoVRadius, self.FoVRadius),
                                      frame = self.frame,
                                      axes  = [self.AxisEnergyReco]
                                      )
            
        elif self.analysis=="1D":
            
            # Center 3D Geometry on Target Coordinates
            # Define Spatial Dimension by Aperture Photometry Region Size
            
            self.RegionRadius = float(configuration['Geometry']['Space']['RegionRadius']) * self.frameUnit
            
            geometry = RegionGeom.create(CircleSkyRegion(self.target, self.RegionRadius),
                                         axes = [self.AxisEnergyReco]
                                         )
        self.geometry = geometry
        
        return None
