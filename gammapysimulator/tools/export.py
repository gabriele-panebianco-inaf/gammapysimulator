#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import os
import pathlib
import shutil

from gammapy.datasets import Datasets
from gammapy.utils.table import table_from_row_data
from matplotlib import use

from gammapysimulator.configure.configure import SimulationConfigurator
from gammapysimulator.tools import utils


class ExportSimulations:
    """
    Class that exports simualation results.
    """
    
    def __init__(self, conf : SimulationConfigurator, datasets, backend='Agg', plotformat='png') -> None:
        """
        Instantiate export class by setting configurator and datasets.
        
        Parameters
        ----------
        conf : gammapysimulator.configure.configure.SimulationConfigurator
            Configurator object.
        datasets : gammapy.datasets.Datasets()
            Simulated datasets.
        backend : str
            Backend to save plots.
        plotformat : str
            Format to save plots.
        """
        # Set attributes
        self.conf= conf
        self.log = conf.log
        self.datasets = datasets
        
        # Add time start and stop arrays
        time_start = [d.gti.time_start.value[0] - d.gti.time_ref.value for d in self.datasets] * u.d
        time_stop  = [d.gti.time_stop.value[0]  - d.gti.time_ref.value for d in self.datasets] * u.d
        
        self.time_start = time_start.to('s')
        self.time_stop  = time_stop.to( 's')
        
        # Set Graphical options
        use(backend)
        self.plotformat=plotformat
        
        return None
    
    def WriteResults(self):
        """
        Write Results according to requested analysis.
        """
        # Make extra results directories
        os.makedirs(self.conf.OutputDirectory.joinpath("plots"))
        os.makedirs(self.conf.OutputDirectory.joinpath("datasets"))
        
        if self.conf.product=="DL4":
            if self.conf.analysis=="1D":
                self.WriteDL4InfoTable()
                self.WriteDL4InfoTable(cumulative=True)
                self.PlotLightCurve()
                self.WriteDatasets()
                self.WriteStacked()
                self.PlotSpectrum()
        elif self.conf.product=="DL3":
            raise NotImplementedError
    
    def WriteDatasets(self):
        """
        Write all datasets.
        """
        self.log.info(f"Write datasets...")
        self.datasets.write(self.conf.OutputDirectory.joinpath("datasets/datasets.fits"),
                            filename_models=self.conf.OutputDirectory.joinpath("models.yaml")
                            )

        return None
    
    def WriteDL4InfoTable(self, cumulative=False):
        """
        Write info table for the datasets.
        
        Parameters
        ----------
        cumulative : bool
            Write Differential or cumulative results
        """
        # Get table
        #info_table = self.datasets.info_table()
        info_table = utils.info_table(self.datasets)
        
        # Add time start and stop columns
        info_table['time_start'] = self.time_start
        info_table['time_stop' ] = self.time_stop
        
        # Write light curve table
        if cumulative:
            table_name = self.conf.OutputDirectory.joinpath('datasets_cumulative.ecsv')
        else:
            table_name = self.conf.OutputDirectory.joinpath('datasets.ecsv')
        
        self.log.info(f"Write {table_name}")
        info_table.write(table_name)
        
        return None
    
    def PlotLightCurve(self, cumulative=False):
        """
        Plot Counts Light Curve.
        
        Parameters
        ----------
        cumulative : bool
            Plot Differential or cumulative results.
        """
        #info_table = self.datasets.info_table(cumulative=cumulative)
        info_table = utils.info_table(self.datasets, cumulative=cumulative)
        
        # Define Quantitites
        counts = info_table['counts']
        counts_errors = np.sqrt(counts)
        background = info_table['background']
        background_errors = np.sqrt(background)
        excess = info_table['excess']
        excess_errors = np.sqrt(np.power(counts_errors,2)+np.power(background_errors,2) )
        time_center =(self.time_stop + self.time_start) / 2.0
        time_errors =(self.time_stop - self.time_start) / 2.0
        e_min = self.conf.axis_energy_reco.bounds[0].value
        e_max = self.conf.axis_energy_reco.bounds[1].value
        
        # Make Plot
        fig, ax = plt.subplots(1, figsize = (10, 5), constrained_layout=True)
        
        # Plot Simulated quantitites
        ax.errorbar(time_center.value, counts, xerr=time_errors.value, yerr=counts_errors,
                    fmt = "o", capsize=2, label = f"Counts")
        ax.errorbar(time_center.value, background, xerr=time_errors.value, yerr=background_errors,
                    fmt = "o", capsize=2, label = f"Background")
        ax.errorbar(time_center.value, excess, xerr=time_errors.value, yerr=excess_errors,
                    fmt = "o", capsize=2, label = f"Excess")
        
        # Plot theoretical quantities
        ax.plot(time_center, info_table['npred'], label="Npred")
        ax.plot(time_center, info_table['npred_signal'], label="Npred Signal")
        ax.plot(time_center, info_table['npred_background'], label="Npred Background")
        
        ax.set_xlabel('Time / s', fontsize = 'large')
        ax.set_ylabel('Counts', fontsize = 'large')
        
        title=f"Counts light curve in [{e_min:2f}, {e_max:2f}] {self.conf.energyUnit}. T0={self.conf.timeRef}."
        ax.set_title(title, fontsize = 'large')
        ax.grid()
        ax.legend(bbox_to_anchor=(1.2, 1.0), loc="upper right")
        
        # Save Plot
        figure_name = self.conf.OutputDirectory.joinpath(f"plots/counts.{self.plotformat}")
        self.log.info(f"Write {figure_name}")
        fig.savefig(figure_name)
        return None
    
    def WriteStacked(self):
        """
        Write Stacked Dataset to get the simulated spectrum.
        """
        stacked = Datasets(self.datasets).stack_reduce(name="stacked")
        stacked.models = self.datasets[0].models
        
        # Write FITS file
        stacked_name = self.conf.OutputDirectory.joinpath(f"datasets/stacked.fits")
        self.log.info(f"Write {stacked_name}")
        stacked.write(stacked_name)
        
        # Write ECSV recap
        stacked_dict = stacked.info_dict()
        stacked_table = table_from_row_data(rows=[stacked_dict])
        stacked_table['time_start'] = self.time_start[0]
        stacked_table['time_stop' ] = self.time_stop[-1]
        table_name = self.conf.OutputDirectory.joinpath(f"stacked.ecsv")
        self.log.info(f"Write {table_name}")
        stacked_table.write(table_name)
        
        return None
    
    def PlotSpectrum(self):
        """
        Plot Spectrum of the Stacked Dataset.
        """
        
        stacked = Datasets(self.datasets).stack_reduce(name="stacked")
        stacked.models = self.datasets[0].models

        # Define Quantitites
        counts = np.squeeze(stacked.counts.data)
        counts_errors = np.sqrt(counts)
        background = np.squeeze(stacked.background.data)
        background_errors = np.sqrt(background)
        excess = np.squeeze(stacked.excess.data)
        excess_errors = np.sqrt(np.power(counts_errors,2)+np.power(background_errors,2) )
        energy_center = stacked.geoms['geom'].axes['energy'].center.value
        energy_errors = stacked.geoms['geom'].axes['energy'].bin_width.value
        npred_counts = np.squeeze(stacked.npred().data)
        npred_signal = np.squeeze(stacked.npred_signal().data)
        npred_background = np.squeeze(stacked.npred_background().data)
        
        time_min = self.time_start[0].value
        time_max = self.time_stop[-1].value
        
        # Make Plot
        fig, ax = plt.subplots(1, figsize = (10, 5), constrained_layout=True)
        
        # Plot Simulated quantitites
        ax.errorbar(energy_center, counts, xerr=energy_errors, yerr=counts_errors,
                    fmt = "o", capsize=2, label = f"Counts")
        ax.errorbar(energy_center, background, xerr=energy_errors, yerr=background_errors,
                    fmt = "o", capsize=2, label = f"Background")
        ax.errorbar(energy_center, excess, xerr=energy_errors, yerr=excess_errors,
                    fmt = "o", capsize=2, label = f"Excess")
        
        # Plot theoretical quantities
        ax.plot(energy_center, npred_counts, label="Npred")
        ax.plot(energy_center, npred_signal, label="Npred Signal")
        ax.plot(energy_center, npred_background, label="Npred Background")

        ax.set_xlabel(f"Energy / {stacked.geoms['geom'].axes['energy'].unit}", fontsize = 'large')
        ax.set_ylabel('Counts', fontsize = 'large')
        ax.set_xscale('log')
        ax.set_yscale('log')
        
        title=f"Counts energy distribution in [{time_min:2f}, {time_max:2f}] {u.s}. T0={self.conf.timeRef}."
        ax.set_title(title, fontsize = 'large')
        ax.grid()
        ax.legend(bbox_to_anchor=(1.2, 1.0), loc="upper right")
        
        # Save Plot
        figure_name = self.conf.OutputDirectory.joinpath(f"plots/spectrum.{self.plotformat}")
        self.log.info(f"Write {figure_name}")
        fig.savefig(figure_name)
        return None
