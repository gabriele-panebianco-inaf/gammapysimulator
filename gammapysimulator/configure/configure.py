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

import yaml

from gammapysimulator.tools import logger, utils

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
        self.log.info("Loading YAML file...")
        with open(ConfigurationFileName) as f:
            configuration = yaml.load(f, Loader = yaml.SafeLoader)
        
        # Simulation Parameters
        self.seed = int(configuration['Simulation']['RandomSeed'])
        
        self.simN = int(configuration['Simulation']['SimuNumber'])
        
        self.OutputDirectory = Path(utils.get_absolute_path(configuration['Simulation']['OutputDirectory'])).absolute()
        self.OutputDirectory.mkdir(parents=True, exist_ok=True)
        
        self.OutputID = str(configuration['Simulation']['OutputID'])
        
        self.product = str(configuration['Simulation']['Product'])
        if self.product not in ["DL3", "DL4"]:
            raise ValueError(f"Only DL3 or DL4 allowed, not {self.product}")
        
        self.analysis = str(configuration['Simulation']['Analysis'])
        if self.analysis not in ["1D", "3D"]:
            raise ValueError(f"Only 1D or 3D allowed, not {self.analysis}")
        
        # Model filename
        self.modelfilename = Path(utils.get_absolute_path(configuration['Model']['FilePath'])).absolute()
        if not self.modelfilename.is_file():
            raise ValueError(f"Model file not found: {self.modelfilename}")
        
        # IRF Parameters
        self.instrument = str(configuration['IRF']['Instrument'])
        if self.instrument not in ["CTA", "Fermi-GBM"]:
            raise ValueError(f"Only CTA or Fermi-GBM allowed, not {self.instrument}")
        
        self.IRFfilename = Path(utils.get_absolute_path(configuration['IRF']['FilePath'])).absolute()
        if not self.IRFfilename.is_file():
            raise ValueError(f"IRF file not found: {self.IRFfilename}")
        
        # Geometry - Time Parameters
        self.timeUnit = u.Unit(str(configuration['Geometry']['Time']['Unit']))
        self.timeStart= float(configuration['Geometry']['Time']['Start'])* self.timeUnit
        self.timeStop = float(configuration['Geometry']['Time']['Stop']) * self.timeUnit
        self.timeReso = float(configuration['Geometry']['Time']['Resolution'])* self.timeUnit
        self.timeRef  = Time(str(configuration['Geometry']['Time']['Reference']))
        
        # Geometry - Energy Axes
        self.energyUnit = u.Unit(str(configuration['Geometry']['Energy']['Unit']))
        
        self.axis_energy_reco = MapAxis.from_energy_bounds(
            configuration['Geometry']['Energy']['RangeReco'][0] * self.energyUnit,
            configuration['Geometry']['Energy']['RangeReco'][1] * self.energyUnit,
            configuration['Geometry']['Energy']['RecoBinPerDecade'],
            per_decade=True,
            name='energy'
            )
        
        self.axis_energy_true = MapAxis.from_energy_bounds(
            configuration['Geometry']['Energy']['RangeTrue'][0] * self.energyUnit,
            configuration['Geometry']['Energy']['RangeTrue'][1] * self.energyUnit,
            configuration['Geometry']['Energy']['TrueBinPerDecade'],
            per_decade=True,
            name='energy_true'
            )
        
        # Geometry - Space Parameters
        self.frame = str(configuration['Geometry']['Space']['Frame'])
        self.frameUnit = u.Unit(str(configuration['Geometry']['Space']['Unit']))
        self.FoVRadius = float(configuration['Geometry']['Space']['FoVRadius']) * self.frameUnit
        self.resolution= float(configuration['Geometry']['Space']['Resolution'])* self.frameUnit
        
        self.pointing = SkyCoord(float(configuration['Geometry']['Space']['PointingLon']),
                                 float(configuration['Geometry']['Space']['PointingLat']),
                                 unit = self.frameUnit,
                                 frame= self.frame
                                 )
        
        if self.analysis=="3D":
            geometry = WcsGeom.create(skydir= self.pointing,
                                      binsz = self.resolution.to('deg').value,
                                      width = (self.FoVRadius, self.FoVRadius),
                                      frame = self.frame,
                                      axes  = [self.axis_energy_reco]
                                      )
        elif self.analysis=="1D":
            geometry = RegionGeom.create(CircleSkyRegion(self.pointing, self.FoVRadius),
                                         axes = [self.axis_energy_reco]
                                         )
        self.geometry = geometry
        
        return None