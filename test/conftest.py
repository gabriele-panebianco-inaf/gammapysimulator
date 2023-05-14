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
from pathlib import Path

@pytest.fixture
def path_repository():
    """Define the absolute path of the repository."""
    return Path(__file__).absolute().parents[1]

@pytest.fixture
def path_configuration_files(path_repository):
    """
    Define the absolute paths of the configuration files.
    """
    
    path_data = path_repository.joinpath("test/data")
    files={'configuration':str(path_data.joinpath("configuration.yml")),
           'model':str(path_data.joinpath("model.yml")),
           'irf':str(path_data.joinpath('Prod5-North-20deg-AverageAz-4LSTs09MSTs.1800s-v0.1.fits'))
           }
    return files

@pytest.fixture
def Mock_Observations(path_configuration_files):
    """Define Mock Observations"""
    
    ObservationsStart= np.linspace(0,4,4,endpoint=False)
    ObservationsStop = ObservationsStart+1.0
    ObservationsStart= ObservationsStart.tolist()* u.h
    ObservationsStop = ObservationsStop.tolist() * u.h
    
    pointing = SkyCoord(83.63, 22.41, frame="fk5", unit="deg")
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
