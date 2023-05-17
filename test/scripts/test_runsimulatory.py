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
    
    def test_Simulate(self, path_repository, path_configuration_files):
        """Test that the class is correctly instantiated."""
        
        os.system(f"runsimulator.py -conf {path_configuration_files['configurationCTA1D']} ")
        OutputDirectory = path_repository.joinpath("SIMULATIONS/TestCTA1D")
        
        # Assert file existence
        assert os.path.isfile(OutputDirectory.joinpath("simulator.log"))
        assert os.path.isfile(OutputDirectory.joinpath("lightcurve.ecsv"))
        assert os.path.isfile(OutputDirectory.joinpath("lightcurve_cumulative.ecsv"))
        assert os.path.isfile(OutputDirectory.joinpath("models.yaml"))
        assert os.path.isfile(OutputDirectory.joinpath("models_covariance.dat"))
        assert os.path.isfile(OutputDirectory.joinpath("stacked_counts.ecsv"))
        assert os.path.isfile(OutputDirectory.joinpath("stacked_spectrum.ecsv"))
        assert os.path.isfile(OutputDirectory.joinpath("plots/lightcurve.png"))
        assert os.path.isfile(OutputDirectory.joinpath("plots/stacked_spectrum.png"))
        assert os.path.isfile(OutputDirectory.joinpath("datasets/datasets.fits"))
        
        for name in ["stacked", "pha_obsonoff-0", "pha_obsonoff-1","pha_obsonoff-2", "pha_obsonoff-3"]:
            assert os.path.isfile(OutputDirectory.joinpath(f"datasets/{name}.fits"))
            assert os.path.isfile(OutputDirectory.joinpath(f"datasets/{name}_arf.fits"))
            assert os.path.isfile(OutputDirectory.joinpath(f"datasets/{name}_rmf.fits"))
            assert os.path.isfile(OutputDirectory.joinpath(f"datasets/{name}_bkg.fits"))
        
        return None