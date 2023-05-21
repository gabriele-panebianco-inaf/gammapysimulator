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
    
    def test_SimulateCTA1D(self, path_repository, path_configuration_files):
        """Test that the script produces the expected files for a CTA 1D Simulation."""
        
        os.system(f"runsimulator.py -conf {path_configuration_files['configurationCTA1D']} ")
        OutputDirectory = path_repository.joinpath("SIMULATIONS/TestCTA1D")
        
        # Assert file existence
        assert os.path.isfile(OutputDirectory.joinpath("simulator.log"))
        assert os.path.isfile(OutputDirectory.joinpath("configuration.yaml"))
        assert os.path.isfile(OutputDirectory.joinpath("models.yaml"))
        assert os.path.isfile(OutputDirectory.joinpath("models_covariance.dat"))
        
        assert os.path.exists(OutputDirectory.joinpath("quicklook/lightcurve.ecsv"))
        assert os.path.exists(OutputDirectory.joinpath("quicklook/lightcurve_cumulative.ecsv"))
        assert os.path.exists(OutputDirectory.joinpath("quicklook/stacked_counts.ecsv"))
        assert os.path.exists(OutputDirectory.joinpath("quicklook/stacked_spectrum.ecsv"))
        assert os.path.isfile(OutputDirectory.joinpath("quicklook/lightcurve.png"))
        assert os.path.isfile(OutputDirectory.joinpath("quicklook/stacked_spectrum.png"))
        
        assert os.path.isfile(OutputDirectory.joinpath("irfs/effectivearea.png"))
        #assert os.path.isfile(OutputDirectory.joinpath("irfs/background.png"))
        #assert os.path.isfile(OutputDirectory.joinpath("irfs/energydispersion.png"))
        #assert os.path.isfile(OutputDirectory.joinpath("irfs/psf.png"))
        
        assert os.path.isfile(OutputDirectory.joinpath("datasets/datasets.fits"))
        assert os.path.isfile(OutputDirectory.joinpath(f"datasets/stacked.fits"))
        #assert os.path.isfile(OutputDirectory.joinpath(f"datasets/stacked_arf.fits"))
        #assert os.path.isfile(OutputDirectory.joinpath(f"datasets/stacked_rmf.fits"))
        #assert os.path.isfile(OutputDirectory.joinpath(f"datasets/stacked_bkg.fits"))
        for i in range(4):
            assert os.path.isfile(OutputDirectory.joinpath(f"datasets/dataset-{i}.fits"))
            #assert os.path.isfile(OutputDirectory.joinpath(f"datasets/pha_obsonoff-{i}_arf.fits"))
            #assert os.path.isfile(OutputDirectory.joinpath(f"datasets/pha_obsonoff-{i}_rmf.fits"))
            #assert os.path.isfile(OutputDirectory.joinpath(f"datasets/pha_obsonoff-{i}_bkg.fits"))
        
        return None
    
    def test_SimulateGBM1D(self, path_repository, path_configuration_files):
        """Test that the script produces the expected files for a GBM 1D Simulation."""
        
        os.system(f"runsimulator.py -conf {path_configuration_files['configurationGBM1D']} ")
        OutputDirectory = path_repository.joinpath("SIMULATIONS/TestGBM1D")
        
        # Assert file existence
        assert os.path.isfile(OutputDirectory.joinpath("simulator.log"))
        assert os.path.isfile(OutputDirectory.joinpath("configuration.yaml"))
        assert os.path.isfile(OutputDirectory.joinpath("models.yaml"))
        assert os.path.isfile(OutputDirectory.joinpath("models_covariance.dat"))
        
        assert os.path.exists(OutputDirectory.joinpath("quicklook/lightcurve.ecsv"))
        assert os.path.exists(OutputDirectory.joinpath("quicklook/lightcurve_cumulative.ecsv"))
        assert os.path.exists(OutputDirectory.joinpath("quicklook/stacked_counts.ecsv"))
        assert os.path.exists(OutputDirectory.joinpath("quicklook/stacked_spectrum.ecsv"))
        assert os.path.isfile(OutputDirectory.joinpath("quicklook/lightcurve.png"))
        assert os.path.isfile(OutputDirectory.joinpath("quicklook/stacked_spectrum.png"))
        
        assert os.path.isfile(OutputDirectory.joinpath("irfs/BackgroundCountsSpectrum.png"))
        assert os.path.isfile(OutputDirectory.joinpath("irfs/BackgroundRates_bak.png"))
        assert os.path.isfile(OutputDirectory.joinpath("irfs/DetectorResponseMatrix_rsp.png"))
        assert os.path.isfile(OutputDirectory.joinpath("irfs/EffectiveArea_rsp.png"))
        assert os.path.isfile(OutputDirectory.joinpath("irfs/EnergyDispersionProbability_rsp.png"))
        assert os.path.isfile(OutputDirectory.joinpath("irfs/Exposure.png"))
        
        assert os.path.isfile(OutputDirectory.joinpath(f"datasets/datasets.fits"))
        assert os.path.isfile(OutputDirectory.joinpath(f"datasets/stacked.fits"))
        #assert os.path.isfile(OutputDirectory.joinpath(f"datasets/stacked_arf.fits"))
        #assert os.path.isfile(OutputDirectory.joinpath(f"datasets/stacked_rmf.fits"))
        #assert os.path.isfile(OutputDirectory.joinpath(f"datasets/stacked_bkg.fits"))
        for i in range(25):
            assert os.path.isfile(OutputDirectory.joinpath(f"datasets/dataset-{i}.fits"))
            #assert os.path.isfile(OutputDirectory.joinpath(f"datasets/pha_obsonoff-{i}_arf.fits"))
            #assert os.path.isfile(OutputDirectory.joinpath(f"datasets/pha_obsonoff-{i}_rmf.fits"))
            #assert os.path.isfile(OutputDirectory.joinpath(f"datasets/pha_obsonoff-{i}_bkg.fits"))
        
        return None
    
    #TODO: Add test for CTA 3D