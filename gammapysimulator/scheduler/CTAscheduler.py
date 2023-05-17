#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

from gammapy.data import Observation, Observations
from gammapy.irf import load_cta_irfs
from gammapysimulator.scheduler.scheduler import Scheduler

import numpy as np

class CTAScheduler(Scheduler):
    """
    Class that organises Observations and Datasets with CTA IRFs.
    """
    
    def __init__(self, configurator) -> None:
        """
        Instantiate Scheduler by reading IRFs.
        
        Parameters
        ----------
        configurator : gammapysimulator.configure.configure.SimulationConfigurator
            Configurator object.
        """
        
        self.conf= configurator
        self.log = configurator.log
        
        # Load IRFs
        self.log.info(f"Load CTA IRFs from file: {self.conf.IRFfilepath}")      
        self.irfs = load_cta_irfs(self.conf.IRFfilepath)


    def SetObservations(self):
        """
        Define the Observations object.
        """
        
        # Create Schedule
        ObservationsStart, ObservationsStop = self.DefineSchedule()
        
        # Define observations
        observations = Observations()
        
        # Observations are time bins, they differ by ID, start and stop
        for i, (tstart, tstop) in enumerate(zip(ObservationsStart, ObservationsStop)):
            obs = Observation.create(self.conf.pointing,
                                     obs_id=i,
                                     tstart=tstart,
                                     tstop=tstop,
                                     irfs=self.irfs,
                                     reference_time=self.conf.timeRef
                                     )
            observations.append(obs)
            
        self.observations = observations
        