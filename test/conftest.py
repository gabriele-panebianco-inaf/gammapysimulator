#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

import pytest

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
