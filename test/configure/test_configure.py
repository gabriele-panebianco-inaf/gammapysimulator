#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

import pytest
from gammapysimulator.configure.configure import SimulationConfigurator

#######################################################
# Fixtures


#######################################################
# Test

@pytest.mark.configure
class TestSimulationConfigurator:
    """Class to test Simulation Configurator"""
    
    def test_init(self):
        """Test that the class is corectly instantiated"""
        configurator = SimulationConfigurator("filename")
        
        assert isinstance(configurator, SimulationConfigurator)
