#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

from gammapysimulator.tools import logger
import yaml
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
        
        return