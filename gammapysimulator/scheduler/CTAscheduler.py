#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

from gammapy.data import Observation, Observations
from gammapy.irf import load_cta_irfs

import numpy as np

class CTAScheduler:
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



    def DefineSchedule(self):
        """
        Define the Schedule: start and stop of each Observation (tme bin).
        """

        # Estimate Number of observations to use numpy linspace
        ObservationsNumber = (self.conf.timeStop - self.conf.timeStart) / self.conf.timeReso
        ObservationsNumber = int(np.floor(ObservationsNumber.to('').value))

        # Define starting time of each observation linearly spaced during the night (wrt reference time)
        ObservationsStart, Livetimes = np.linspace(self.conf.timeStart.value,
                                                   self.conf.timeStop.value,
                                                   num = ObservationsNumber,
                                                   endpoint=False,
                                                   retstep=True
                                                   )
        # Compute Stops and Number
        ObservationsStop = ObservationsStart+Livetimes
        ObservationsNumber = ObservationsStart.size
        
        # Convert to Astropy Quantity
        ObservationsStart= ObservationsStart.tolist()* self.conf.timeUnit
        ObservationsStop = ObservationsStop.tolist() * self.conf.timeUnit
        Livetimes = Livetimes * self.conf.timeUnit
        
        # Log
        self.log.info(f"Number of Observations: {ObservationsNumber}. Livetimes: {Livetimes}.")
        self.log.info(f"Observation interval: [{ObservationsStart[0]},{ObservationsStop[-1]}].")

        return ObservationsStart, ObservationsStop
    
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
        