#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

from abc import ABC, abstractmethod
import numpy as np

class Scheduler(ABC):
    """
    Abstract class for different instrument schedulers.
    """
    
    def DefineSchedule(self):
        """
        Define the Schedule: start and stop of each Observation (time bin).
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
    
    def getEmptyDatasets(self):
        """
        Return attribute emptydatasets that contains Datasets with GTI and Reduced IRFs.
        """
        return self.emptydatasets
    
    @abstractmethod
    def ScheduleDatasets(self):
        pass

    