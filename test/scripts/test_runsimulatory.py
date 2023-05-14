#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

import pytest
import os

from gammapysimulator.scripts.runsimulator import Simulate

#######################################################
# Fixtures

#######################################################
# Test

@pytest.mark.scripts
class TestRunSimulator:
    """Class to test RunSimulator"""
    
    def test_Simulate(self, path_configuration_files):
        """Test that the class is correctly instantiated."""
        
        os.system(f"runsimulator.py -conf {path_configuration_files['configuration']} ")
        
        #TODO: assert products