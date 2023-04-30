#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

import pytest

from gammapysimulator.configure import configure
from gammapysimulator.simulator import simulator
from gammapysimulator.tools import logger

#######################################################
# Fixtures


#######################################################
# Test

@pytest.mark.simulator
class TestSimulator:
    """Class to test Simulator"""
    
    def test_init(self, path_configuration_files):
        """Test that the class is correctly instantiated"""
        
        # Instantiate Configurator
        sourcesimulator = simulator.Simulator(path_configuration_files['configuration'])
        
        # Assert Class
        assert isinstance(sourcesimulator, simulator.Simulator)
        assert isinstance(sourcesimulator.conf, configure.SimulationConfigurator)
        