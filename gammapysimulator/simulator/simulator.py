#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

from gammapysimulator.configure.configure import SimulationConfigurator


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
        
        # Read Configuration file
        configurator = SimulationConfigurator()
        configurator.read(ConfigurationFileName)
        
        self.conf = configurator