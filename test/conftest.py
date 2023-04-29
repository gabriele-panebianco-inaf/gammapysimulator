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
    files={'configuration':str(path_data.joinpath("configuration.yml"))
           }
    return files
