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
from gammapy.maps import MapAxis, WcsGeom
from pathlib import Path

from gammapysimulator.configure import configure
from gammapysimulator.scheduler import CTAscheduler

#######################################################
# Fixtures
@pytest.fixture
def mock_Observations(path_configuration_files):
    """Define Mock Observations"""
    
    ObservationsStart= np.linspace(-5,20,125,endpoint=False)
    ObservationsStop = ObservationsStart+0.2
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

@pytest.mark.scheduler
class TestCTAScheduler:
    """Class to test CTA Scheduler"""
    
    def test_init(self, path_configuration_files):
        """Test Scheduler initialization"""
        
        # Prepare Configurator
        configurator = configure.SimulationConfigurator()
        configurator.read(path_configuration_files['configuration'])
        
        # Define Scheduler
        scheduler = CTAscheduler.CTAScheduler(configurator)
        
        # Assert Class
        assert isinstance(scheduler, CTAscheduler.CTAScheduler)
        assert isinstance(scheduler.conf, configure.SimulationConfigurator)
        assert "aeff"  in scheduler.irfs
        assert "psf"   in scheduler.irfs
        assert "edisp" in scheduler.irfs
        assert "bkg"   in scheduler.irfs
        
    def test_DefineSchedule(self, path_configuration_files):
        """Test function to compute observation start and stop times"""
        
        # Define scheduler
        configurator = configure.SimulationConfigurator()
        configurator.read(path_configuration_files['configuration'])
        scheduler = CTAscheduler.CTAScheduler(configurator)
        
        # Create Schedule
        ObservationsStart, ObservationsStop = scheduler.DefineSchedule()
        
        # Assert size
        assert ObservationsStart.size==125
        assert ObservationsStop.size==125
        
        # Assert content
        my_start = np.linspace(-5,20,125,endpoint=False)
        my_stop = my_start+0.2
        assert not np.any(ObservationsStart.value-my_start)
        assert not np.any(ObservationsStop.value-my_stop)
        
    def test_SetObservations(self, path_configuration_files, mock_Observations):
        """Test Setting Observations"""
        
        # Define scheduler
        configurator = configure.SimulationConfigurator()
        configurator.read(path_configuration_files['configuration'])
        scheduler = CTAscheduler.CTAScheduler(configurator)
        scheduler.SetObservations()
        
        # Test Number of Observations
        assert len(scheduler.observations)==len(mock_Observations)
        
        # Test Observations attributes
        for sched_obs, mock_obs in zip(scheduler.observations, mock_Observations):
            assert sched_obs.aeff == mock_obs.aeff
            assert sched_obs.available_hdus == mock_obs.available_hdus
            assert sched_obs.available_irfs == mock_obs.available_irfs
            assert sched_obs.bkg == mock_obs.bkg
            assert sched_obs.edisp == mock_obs.edisp
            assert sched_obs.muoneff == mock_obs.muoneff
            assert sched_obs.obs_info == mock_obs.obs_info
            assert sched_obs.observation_dead_time_fraction == mock_obs.observation_dead_time_fraction
            assert sched_obs.observation_time_duration == mock_obs.observation_time_duration
            assert sched_obs.pointing_radec == mock_obs.pointing_radec
            assert sched_obs.psf == mock_obs.psf
            assert sched_obs.rad_max == mock_obs.rad_max
            assert sched_obs.tstart == mock_obs.tstart
            assert sched_obs.tstop == mock_obs.tstop
            assert sched_obs.gti.time_delta == mock_obs.gti.time_delta
            assert sched_obs.gti.time_intervals == mock_obs.gti.time_intervals
            assert sched_obs.gti.time_ref == mock_obs.gti.time_ref
            assert sched_obs.gti.time_start == mock_obs.gti.time_start
            assert sched_obs.gti.time_stop == mock_obs.gti.time_stop
            assert sched_obs.gti.time_sum == mock_obs.gti.time_sum

