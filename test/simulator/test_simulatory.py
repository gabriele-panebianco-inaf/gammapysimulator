#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

import numpy as np
import os
import pytest

from gammapy.data import Observations
from gammapy.datasets import Datasets
from gammapy.modeling.models import SkyModel, Models

from gammapysimulator.configure import configure
from gammapysimulator.simulator import simulator
from gammapysimulator.tools import utils

#######################################################
# Fixtures

#######################################################
# Test

@pytest.mark.simulator
class TestSimulatorCTA:
    """Class to test Simulator with CTA IRFs."""
    
    def test_init(self, path_configuration_files):
        """Test that the class is correctly instantiated."""
        
        # Instantiate Simulator
        sourcesimulator = simulator.Simulator(path_configuration_files['configurationCTA1D'])
        
        # Assert Class and attributes
        assert isinstance(sourcesimulator, simulator.Simulator)
        assert isinstance(sourcesimulator.conf, configure.SimulationConfigurator)
        assert isinstance(sourcesimulator.models, Models)
        assert isinstance(sourcesimulator.observations, Observations)
        assert sourcesimulator.psf_containment
        
    def test_SetObservations(self, path_configuration_files, Mock_Observations):
        """Test that Observations are correctly read."""
        
        # Instantiate Configurator
        sourcesimulator = simulator.Simulator(path_configuration_files['configurationCTA1D'])
        observations = sourcesimulator.SetObservations()
        
        # Assert PSF Containment
        assert sourcesimulator.psf_containment
        
        # Test Number of Observations
        assert len(sourcesimulator.observations)==len(Mock_Observations)
        
        # Test Observations attributes
        for simu_obs, mock_obs in zip(sourcesimulator.observations, Mock_Observations):
            assert simu_obs.aeff == mock_obs.aeff
            assert simu_obs.available_hdus == mock_obs.available_hdus
            assert simu_obs.available_irfs == mock_obs.available_irfs
            assert simu_obs.bkg == mock_obs.bkg
            assert simu_obs.edisp == mock_obs.edisp
            assert simu_obs.muoneff == mock_obs.muoneff
            assert simu_obs.obs_info == mock_obs.obs_info
            assert simu_obs.observation_dead_time_fraction == mock_obs.observation_dead_time_fraction
            assert simu_obs.observation_time_duration == mock_obs.observation_time_duration
            assert simu_obs.pointing_radec == mock_obs.pointing_radec
            assert simu_obs.psf == mock_obs.psf
            assert simu_obs.rad_max == mock_obs.rad_max
            assert simu_obs.tstart == mock_obs.tstart
            assert simu_obs.tstop == mock_obs.tstop
            assert simu_obs.gti.time_delta == mock_obs.gti.time_delta
            assert simu_obs.gti.time_intervals == mock_obs.gti.time_intervals
            assert simu_obs.gti.time_ref == mock_obs.gti.time_ref
            assert simu_obs.gti.time_start == mock_obs.gti.time_start
            assert simu_obs.gti.time_stop == mock_obs.gti.time_stop
            assert simu_obs.gti.time_sum == mock_obs.gti.time_sum
            
    def test_Simulate1D(self, path_configuration_files):
        """Test that Datasets are correctly set."""
        
        # Instantiate Simulator
        sourcesimulator = simulator.Simulator(path_configuration_files['configurationCTA1D'])
        
        # Run Dataset Reduction
        datasets = sourcesimulator.SimulateCTA()
        
        assert isinstance(datasets, Datasets)
        
        # Test Number of Observations
        assert len(sourcesimulator.observations)==4
        
        # Assert other info about datasets
        assert datasets.energy_axes_are_aligned
        assert datasets.energy_ranges[0][0].value == pytest.approx(0.3)
        assert datasets.energy_ranges[1][0].value == pytest.approx(100)
        
        # Assert information on light curve
        table_diff = utils.info_table(datasets)
        
        assert table_diff['npred_signal'].data == pytest.approx(np.array([407.0,407.0,407.0,407.0]), abs=1.0)
        assert table_diff['npred_background'].data == pytest.approx(np.array([159.0, 159.0, 159.0, 159.0]), abs=7.0)
        assert table_diff['counts'].data == pytest.approx(np.array([566, 583, 534, 574]),abs=1.0)
        assert table_diff['background'].data == pytest.approx(np.array([165.8, 158.8, 161.0, 155.0]),abs=1.0)

        # Assert information on stacked spectrum
        stacked = datasets.stack_reduce(name="stacked")
        stacked.models = datasets[0].models
        
        assert stacked.geoms['geom'].axes['energy'].nbin == 11
        
        assert np.squeeze(stacked.counts.data) == pytest.approx(np.array([810,491,312,257,163,95,73,26,17,9,4]),abs=1.0)
        assert np.squeeze(stacked.background.data)== pytest.approx(np.array([312.2,151.6,87.8,54.4,18.4,8.0,4.4,1.4,1.2,0.4,0.8]),abs=1.0)
        assert np.squeeze(stacked.npred_signal().data)== pytest.approx(np.array([479.7,366.1,254.7,191.4,142.7,93.8,47.9,22.8,14.8,8.9,4.7]),abs=1.0)
        assert np.squeeze(stacked.npred_background().data)== pytest.approx(np.array([313.5,150.1,86.32,54.9,18.4,7.9,4.8,1.4,1.2,0.4,0.8]),abs=1.0)
        
        stacked_dict = stacked.info_dict()
        
        assert stacked_dict['counts'] == 2257
        assert stacked_dict['sqrt_ts'] == pytest.approx(43, abs=1.0)
        assert stacked_dict['background'] == pytest.approx(640, abs=1.0)
        assert stacked_dict['exposure_min'].value == pytest.approx( 2989892864.0, abs=1.0)
        assert stacked_dict['exposure_max'].value == pytest.approx(17710960640.0, abs=1.0)
        assert stacked_dict['livetime'].value == pytest.approx(14400, abs=1.0)
        assert stacked_dict['ontime'].value == pytest.approx(14400, abs=1.0)
        
        return None

    def test_ExportResults1D(self, path_repository, path_configuration_files):
        """Test that expected results are written."""
        
        # Instantiate Simulator
        sourcesimulator = simulator.Simulator(path_configuration_files['configurationCTA1D'])
        simulations = sourcesimulator.RunSimulation()
        sourcesimulator.ExportResults(simulations)
        
        # Assert file existence
        OutputDirectory = sourcesimulator.conf.OutputDirectory
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