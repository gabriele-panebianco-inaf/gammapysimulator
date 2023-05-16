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
    files={'configurationCTA1D':str(path_data.joinpath("configurationCTA1D.yml")),
           'configurationCTA3D':str(path_data.joinpath("configurationCTA3D.yml")),
           'modelCrab':str(path_data.joinpath("modelCrab.yml")),
           'irfCTA':str(path_data.joinpath('Prod5-North-20deg-AverageAz-4LSTs09MSTs.1800s-v0.1.fits'))
           }
    return files

@pytest.fixture
def Mock_Observations(path_configuration_files):
    """Define Mock Observations"""
    
    ObservationsStart= np.linspace(0,4,4,endpoint=False)
    ObservationsStop = ObservationsStart+1.0
    ObservationsStart= ObservationsStart.tolist()* u.h
    ObservationsStop = ObservationsStop.tolist() * u.h
    
    pointing = SkyCoord(83.63, 22.41, frame="icrs", unit="deg")
    irfs = load_cta_irfs(path_configuration_files['irfCTA'])
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
