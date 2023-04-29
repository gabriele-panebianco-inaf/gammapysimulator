#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

from gammapysimulator.tools import logger

class SimulationConfigurator:
    """
    Class that reads a YML file and defines the parameters needed for the simulation.
    """
    
    def __init__(self, ConfigurationFileName)->None:
        
        self.log = logger.SimulatorLogger().getLogger()
        
        self.log.info(f"FileName: {ConfigurationFileName}")
        return None