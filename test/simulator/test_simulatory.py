#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

import pytest

import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from gammapy.data import Observation, Observations
from gammapy.irf import load_cta_irfs
from gammapy.modeling.models import SkyModel

from gammapysimulator.configure import configure
from gammapysimulator.simulator import simulator

#######################################################
# Fixtures
@pytest.fixture
def mock_Observations(path_configuration_files):
    """Define Mock Observations"""
    
    ObservationsStart= np.linspace(-5,20,25,endpoint=False)
    ObservationsStop = ObservationsStart+1.0
    ObservationsStart= ObservationsStart.tolist()* u.s
    ObservationsStop = ObservationsStop.tolist() * u.s
    
    pointing = SkyCoord(83.63, 22.01, frame="fk5", unit="deg")
    irfs = load_cta_irfs(path_configuration_files['irf'])
    timeRef = Time("2023-01-01T00:00:00")
    
    observations = Observations()
    
    for i, (tstart, tstop) in enumerate(zip(ObservationsStart, ObservationsStop)):
        obs = Observation.create(pointing,
                                 obs_id=i,
                                 tstart=tstart,
                                 tstop=tstop,
                                 irfs=irfs,
                                 reference_time=timeRef
                                 )
        observations.append(obs)
        
    return observations


#######################################################
# Test

@pytest.mark.simulator
class TestSimulator:
    """Class to test Simulator"""
    
    def test_init(self, path_configuration_files):
        """Test that the class is correctly instantiated."""
        
        # Instantiate Configurator
        sourcesimulator = simulator.Simulator(path_configuration_files['configuration'])
        
        # Assert Class and attributes
        assert isinstance(sourcesimulator, simulator.Simulator)
        assert isinstance(sourcesimulator.conf, configure.SimulationConfigurator)
        assert isinstance(sourcesimulator.models, SkyModel)
        assert isinstance(sourcesimulator.observations, Observations)
        assert sourcesimulator.psf_containment
        
    def test_SetObservations(self, path_configuration_files, mock_Observations):
        """Test that Observations are correctly read."""
        
        # Instantiate Configurator
        sourcesimulator = simulator.Simulator(path_configuration_files['configuration'])
        observations = sourcesimulator.SetObservations()
        
        # Test Number of Observations
        assert len(sourcesimulator.observations)==len(mock_Observations)
        
        # Test Observations attributes
        for simu_obs, mock_obs in zip(sourcesimulator.observations, mock_Observations):
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
        