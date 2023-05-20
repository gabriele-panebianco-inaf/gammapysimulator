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
        assert isinstance(sourcesimulator.emptydatasets, Datasets)
        assert sourcesimulator.psf_containment
        
    def test_SetDL4IRFs(self, path_configuration_files, Mock_Observations):
        """Test that DL4 Datasets are correctly read."""
        
        # Instantiate Configurator
        sourcesimulator = simulator.Simulator(path_configuration_files['configurationCTA1D'])
        emptydatasets = sourcesimulator.SetDL4IRFs()
        
        # Assert PSF Containment
        assert sourcesimulator.psf_containment
        
        # Test Number of Observations
        assert len(sourcesimulator.emptydatasets)==len(Mock_Observations)
        
        # TODO: Test Datasets attributes
            
    def test_Simulate1D(self, path_configuration_files):
        """Test that Datasets are correctly set."""
        
        # Instantiate Simulator
        sourcesimulator = simulator.Simulator(path_configuration_files['configurationCTA1D'])
        
        # Run Dataset Reduction
        datasets = sourcesimulator.RunSimulation()
        
        assert isinstance(datasets, Datasets)
        
        # Test Number of Datasets
        assert len(datasets)==4
        
        # Assert other info about datasets
        assert datasets.energy_axes_are_aligned
        assert datasets.energy_ranges[0][0].value == pytest.approx(0.3)
        assert datasets.energy_ranges[1][0].value == pytest.approx(100)
        
        # Assert information on light curve
        #table_diff = datasets.info_table()
        table_diff = utils.info_table(datasets)
        
        assert table_diff['npred_signal'].data == pytest.approx(407.0*np.ones(4), abs=1.0)
        assert table_diff['npred_background'].data == pytest.approx(162.4*np.ones(4), abs=1.0)
        assert table_diff['counts'].data == pytest.approx(np.array([534, 574, 548, 578]),abs=1.0)

        # Assert information on stacked spectrum
        stacked = datasets.stack_reduce(name="stacked")
        stacked.models = datasets[0].models
        
        assert stacked.geoms['geom'].axes['energy'].nbin == 11
        
        assert np.squeeze(stacked.counts.data) == pytest.approx(np.array([782, 516, 340, 233, 152, 103, 50, 24, 21, 6, 7]),abs=1.0)
        assert np.squeeze(stacked.npred_signal().data)== pytest.approx(np.array([479.7, 366.1, 254.7, 191.4, 142.7, 93.8, 47.9, 22.8, 14.8, 8.9, 4.7]), abs=1.0)
        assert np.squeeze(stacked.npred_background().data)== pytest.approx(np.array([308.3, 160.8, 93.8, 49.5, 21.0, 8.7, 3.7, 1.5, 1.0, 0.8, 0.6]), abs=1.0)
        
        stacked_dict = stacked.info_dict()
        
        assert stacked_dict['counts'] == 2234
        assert stacked_dict['sqrt_ts'] == pytest.approx(48, abs=1.0)
        assert stacked_dict['background'] == pytest.approx(649, abs=1.0)
        assert stacked_dict['exposure_min'].value == pytest.approx( 2989892864.0, abs=1.0)
        assert stacked_dict['exposure_max'].value == pytest.approx(17710960640.0, abs=1.0)
        assert stacked_dict['livetime'].value == pytest.approx(14400, abs=1.0)
        assert stacked_dict['ontime'].value == pytest.approx(14400, abs=1.0)
        
        return None
