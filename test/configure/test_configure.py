#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

import pytest

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from gammapy.maps import MapAxis, WcsGeom
from pathlib import Path

from gammapysimulator.configure.configure import SimulationConfigurator
from gammapysimulator.tools.logger import SimulatorLogger

#######################################################
# Fixtures


#######################################################
# Test

@pytest.mark.configure
class TestSimulationConfigurator:
    """Class to test Simulation Configurator"""
    
    def test_init(self):
        """Test that the class is correctly instantiated"""
        
        # Instantiate Configurator
        configurator = SimulationConfigurator()
        
        # Assert Class
        assert isinstance(configurator, SimulationConfigurator)
        assert isinstance(configurator.log, SimulatorLogger)

    def test_init(self, path_repository, path_configuration_files):
        """Test that a YAML configuration file is read"""
        
        # Instantiate Configurator
        configurator = SimulationConfigurator()        
        # Read Configuration file
        configurator.read(path_configuration_files['configuration'])
        
        # Assert attributes value
        assert configurator.seed == 7
        assert configurator.simN == 2
        assert configurator.OutputDirectory == path_repository.joinpath("SIMULATIONS/")
        assert configurator.OutputID == "test"
        assert configurator.product == "DL3"
        assert configurator.analysis == "3D"
        assert configurator.modelfilename == Path(path_configuration_files['model'])
        assert configurator.instrument == "CTA"
        assert configurator.IRFfilename== Path(path_configuration_files['irf'])
        assert configurator.timeUnit == u.s
        assert configurator.timeStart.value== pytest.approx(-5,1E-2)
        assert configurator.timeStop.value == pytest.approx(20,1E-2)
        assert configurator.timeReso.value == pytest.approx(1.0,1E-2)
        assert configurator.timeRef  == Time("2023-01-01T00:00:00")
        assert configurator.energyUnit == u.TeV
        assert configurator.axis_energy_reco == MapAxis.from_energy_bounds("0.3 TeV", "100 TeV", 4, per_decade=True,name='energy')
        assert configurator.axis_energy_true == MapAxis.from_energy_bounds("0.1 TeV", "300 TeV", 5, per_decade=True,name='energy_true')
        assert configurator.frame =="fk5"
        assert configurator.frameUnit == u.deg
        assert configurator.FoVRadius.value == pytest.approx(5,1E-2)
        assert configurator.resolution.value== pytest.approx(0.02,1E-2)
        assert configurator.pointing == SkyCoord(83.63, 22.01, frame="fk5", unit="deg")
        assert configurator.geometry == WcsGeom.create(skydir= configurator.pointing,
                                                       binsz = 0.02,
                                                       width = (5*u.deg, 5*u.deg),
                                                       frame = "fk5",
                                                       axes  = [configurator.axis_energy_reco]
                                                       )
        
        
