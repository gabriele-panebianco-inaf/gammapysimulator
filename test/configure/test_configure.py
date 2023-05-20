#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

import os
import pytest

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from gammapy.maps import MapAxis, WcsGeom, RegionGeom
from pathlib import Path
from regions import CircleSkyRegion

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
        #assert isinstance(configurator.log, SimulatorLogger)

    def test_initCTA1D(self, path_repository, path_configuration_files):
        """
        Test that a YAML configuration file is read.
        Use CTA 1D configuration.       
        """
        
        # Instantiate Configurator
        configurator = SimulationConfigurator()        
        # Read Configuration file
        configurator.read(path_configuration_files['configurationCTA1D'])
        
        # Assert log file exists
        assert os.path.isfile(path_repository.joinpath("SIMULATIONS/TestCTA1D/simulator.log"))
        
        # Assert attributes value
        
        # Simulation
        assert configurator.seed == 7
        assert configurator.simN == 2
        assert configurator.OutputDirectory == path_repository.joinpath("SIMULATIONS/TestCTA1D")
        assert configurator.OutputID == "TestCTA1D"
        assert configurator.product == "DL4"
        assert configurator.analysis == "1D"
        
        # Model
        assert configurator.modelfilepath == Path(path_configuration_files['modelCrab'])
        
        # IRF
        assert configurator.instrument == "CTA"
        assert configurator.detector == "North-4LSTs-09MSTs"
        assert configurator.irf_pointlike == False
        assert configurator.IRFfilepath== Path(path_configuration_files['irfCTA'])
        
        # Target and Pointing
        assert configurator.pointing == SkyCoord(83.63, 22.41, frame="icrs", unit="deg")
        assert configurator.target   == SkyCoord(83.63, 22.01, frame="icrs", unit="deg")
        
        # Geometry - Time
        assert configurator.timeUnit == u.h
        assert configurator.timeStart.value== pytest.approx(0.0, 1E-2)
        assert configurator.timeStop.value == pytest.approx(4.0, 1E-2)
        assert configurator.timeReso.value == pytest.approx(1.0, 1E-2)
        assert configurator.timeRef  == Time("2023-01-01T00:00:00")
        
        # Geometry - Energy
        assert configurator.energyUnit == u.TeV
        assert configurator.AxisEnergyReco == MapAxis.from_energy_bounds("0.3 TeV", "100 TeV", 4, per_decade=True,name='energy')
        assert configurator.AxisEnergyTrue == MapAxis.from_energy_bounds("0.1 TeV", "300 TeV", 5, per_decade=True,name='energy_true')
        
        # Geometry - Space
        assert configurator.frame =="icrs"
        assert configurator.frameUnit == u.deg
        assert configurator.RegionRadius.value == pytest.approx(0.2,1E-2)
        assert not hasattr(configurator, 'FoVRadius')
        assert not hasattr(configurator, 'resolution')

        assert configurator.geometry == RegionGeom.create(CircleSkyRegion(configurator.target,
                                                                          configurator.RegionRadius
                                                                          ),
                                                          axes = [configurator.AxisEnergyReco]
                                                          )
        
    def test_initCTA3D(self, path_repository, path_configuration_files):
        """
        Test that a YAML configuration file is read.
        Use CTA 3D configuration.      
        """
        
        # Instantiate Configurator
        configurator = SimulationConfigurator()        
        # Read Configuration file
        configurator.read(path_configuration_files['configurationCTA3D'])
        
        # Assert log file exists
        assert os.path.isfile(path_repository.joinpath("SIMULATIONS/TestCTA3D/simulator.log"))
        
        # Assert attributes value
        
        # Simulation
        assert configurator.seed == 7
        assert configurator.simN == 2
        assert configurator.OutputDirectory == path_repository.joinpath("SIMULATIONS/TestCTA3D")
        assert configurator.OutputID == "TestCTA3D"
        assert configurator.product == "DL4"
        assert configurator.analysis == "3D"
        
        # Model
        assert configurator.modelfilepath == Path(path_configuration_files['modelCrab'])
        
        # IRF
        assert configurator.instrument == "CTA"
        assert configurator.detector == "North-4LSTs-09MSTs"
        assert configurator.irf_pointlike == False
        assert configurator.IRFfilepath== Path(path_configuration_files['irfCTA'])
        
        # Target and Pointing
        assert configurator.pointing == SkyCoord(83.63, 22.41, frame="icrs", unit="deg")
        assert configurator.target   == SkyCoord(83.63, 22.01, frame="icrs", unit="deg")
        
        # Geometry - Time
        assert configurator.timeUnit == u.h
        assert configurator.timeStart.value== pytest.approx(0.0, 1E-2)
        assert configurator.timeStop.value == pytest.approx(4.0, 1E-2)
        assert configurator.timeReso.value == pytest.approx(1.0, 1E-2)
        assert configurator.timeRef  == Time("2023-01-01T00:00:00")
        
        # Geometry - Energy
        assert configurator.energyUnit == u.TeV
        assert configurator.AxisEnergyReco == MapAxis.from_energy_bounds("0.3 TeV", "100 TeV", 4, per_decade=True,name='energy')
        assert configurator.AxisEnergyTrue == MapAxis.from_energy_bounds("0.1 TeV", "300 TeV", 5, per_decade=True,name='energy_true')
        
        # Geometry - Space
        assert configurator.frame =="icrs"
        assert configurator.frameUnit == u.deg
        assert not hasattr(configurator, 'RegionRadius')
        assert configurator.FoVRadius.value == pytest.approx(5.0,1E-2)
        assert configurator.resolution.value == pytest.approx(0.02,1E-2)

        assert configurator.geometry == WcsGeom.create(skydir= configurator.pointing,
                                                       binsz = 0.02,
                                                       width = (5*u.deg, 5*u.deg),
                                                       frame = "icrs",
                                                       axes  = [configurator.AxisEnergyReco]
                                                       )
        
    def test_initGBM1D(self, path_repository, path_configuration_files):
        """
        Test that a YAML configuration file is read.
        Use GBM 1D configuration.       
        """
        
        # Instantiate Configurator
        configurator = SimulationConfigurator()        
        # Read Configuration file
        configurator.read(path_configuration_files['configurationGBM1D'])
        
        # Assert log file exists
        assert os.path.isfile(path_repository.joinpath("SIMULATIONS/TestGBM1D/simulator.log"))
        
        # Assert attributes value
        
        # Simulation
        assert configurator.seed == 7
        assert configurator.simN == 2
        assert configurator.OutputDirectory == path_repository.joinpath("SIMULATIONS/TestGBM1D")
        assert configurator.OutputID == "TestGBM1D"
        assert configurator.product == "DL4"
        assert configurator.analysis == "1D"
        
        # Model
        assert configurator.modelfilepath == Path(path_configuration_files['modelGRB'])
        
        # IRF
        assert configurator.instrument == "Fermi-GBM"
        assert configurator.detector == "n1"
        assert configurator.irf_pointlike == True
        assert configurator.IRFfilepath['RSP'] == Path(path_configuration_files['irfGBMrsp'])
        assert configurator.IRFfilepath['BAK'] == Path(path_configuration_files['irfGBMbak'])
        
        # Target and Pointing
        assert configurator.pointing == SkyCoord(177.52, 40.0, frame="galactic", unit="deg")
        assert configurator.target   == SkyCoord(177.52, 40.0, frame="galactic", unit="deg")
        
        # Geometry - Time
        assert configurator.timeUnit == u.s
        assert configurator.timeStart.value== pytest.approx(-5.0, 1E-2)
        assert configurator.timeStop.value == pytest.approx(20.0, 1E-2)
        assert configurator.timeReso.value == pytest.approx(1.0, 1E-2)
        assert configurator.timeRef  == Time("2016-05-30T16:01:12")
        
        # Geometry - Energy
        assert configurator.energyUnit == u.keV
        assert configurator.AxisEnergyReco == MapAxis.from_energy_bounds( "10 keV",  "900 keV", 10, per_decade=True,name='energy')
        assert configurator.AxisEnergyTrue == MapAxis.from_energy_bounds("5.0 keV", "1000 keV", 20, per_decade=True,name='energy_true')
        
        # Geometry - Space
        assert configurator.frame =="galactic"
        assert configurator.frameUnit == u.deg
        assert configurator.RegionRadius.value == pytest.approx(0.2,1E-2)
        assert not hasattr(configurator, 'FoVRadius')
        assert not hasattr(configurator, 'resolution')

        assert configurator.geometry == RegionGeom.create(CircleSkyRegion(configurator.target,
                                                                          configurator.RegionRadius
                                                                          ),
                                                          axes = [configurator.AxisEnergyReco]
                                                          )
        
        
