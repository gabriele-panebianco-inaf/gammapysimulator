#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

import pytest

import numpy as np

from gammapysimulator.configure import configure
from gammapysimulator.scheduler import CTAscheduler
from gammapysimulator.tools import export

#######################################################
# Fixtures

#######################################################
# Test

@pytest.mark.scheduler
class TestCTAScheduler:
    """Class to test CTA Scheduler"""
    
    def test_init(self, path_configuration_files):
        """Test Scheduler initialization"""
        
        # Prepare Configurator and Export
        configurator = configure.SimulationConfigurator()
        configurator.read(path_configuration_files['configurationCTA1D'])
        exporter = export.ExportSimulations(configurator)
        
        # Define Scheduler
        scheduler = CTAscheduler.CTAScheduler(configurator, exporter)
        
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
        configurator.read(path_configuration_files['configurationCTA1D'])
        exporter = export.ExportSimulations(configurator)
        scheduler = CTAscheduler.CTAScheduler(configurator, exporter)
        
        # Create Schedule
        ObservationsStart, ObservationsStop = scheduler.DefineSchedule()
        
        # Assert size
        assert ObservationsStart.size==4
        assert ObservationsStop.size==4
        
        # Assert content
        my_start = np.linspace(0,4,4,endpoint=False)
        my_stop = my_start+1.0
        assert not np.any(ObservationsStart.value-my_start)
        assert not np.any(ObservationsStop.value-my_stop)
        
    def test_SetObservations(self, path_configuration_files, Mock_Observations):
        """Test Setting Observations"""
        
        # Define scheduler
        configurator = configure.SimulationConfigurator()
        configurator.read(path_configuration_files['configurationCTA1D'])
        exporter = export.ExportSimulations(configurator)
        scheduler = CTAscheduler.CTAScheduler(configurator, exporter)
        observations = scheduler.SetObservations()
        
        # Test Number of Observations
        assert len(observations)==len(Mock_Observations)
        
        # Test Observations attributes
        for sched_obs, mock_obs in zip(observations, Mock_Observations):
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

