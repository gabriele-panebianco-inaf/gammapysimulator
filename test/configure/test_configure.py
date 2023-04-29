#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

import pytest
from gammapysimulator.configure.configure import SimulationConfigurator
from gammapysimulator.tools.logger import SimulatorLogger

#######################################################
# Fixtures


#######################################################
# Test

@pytest.mark.configure
class TestSimulationConfigurator:
    """Class to test Simulation Configurator"""
    
    def test_init(self):
        """Test that the class is correctly instantiated"""
        
        # Instantiate Configurator
        configurator = SimulationConfigurator()
        
        # Assert Class
        assert isinstance(configurator, SimulationConfigurator)
        assert isinstance(configurator.log, SimulatorLogger)

    def test_init(self, path_configuration_files):
        """Test that a YAML configuration file is read"""
        
        # Instantiate Configurator
        configurator = SimulationConfigurator()        
        # Read Configuration file
        configurator.read(path_configuration_files['configuration'])
        
        
        
        
