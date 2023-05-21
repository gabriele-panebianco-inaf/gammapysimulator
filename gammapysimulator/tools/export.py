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

from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table
from gammapy.datasets import Datasets
from gammapy.utils.table import table_from_row_data
from matplotlib import use
from matplotlib.colors import LogNorm, PowerNorm

from gammapysimulator.configure.configure import SimulationConfigurator
from gammapysimulator.tools import utils


class ExportSimulations:
    """
    Class that exports simualation results.
    """
    
    def __init__(self, conf : SimulationConfigurator, backend='Agg', plotformat='png') -> None:
        """
        Instantiate export class by setting configurator and datasets.
        
        Parameters
        ----------
        conf : gammapysimulator.configure.configure.SimulationConfigurator
            Configurator object.
        backend : str
            Backend to save plots.
        plotformat : str
            Format to save plots.
        """
        # Set attributes
        self.conf= conf
        self.log = conf.log
        
        # Make extra directories to store results
        os.makedirs(self.conf.OutputDirectory.joinpath("datasets"))
        os.makedirs(self.conf.OutputDirectory.joinpath("irfs"))
        os.makedirs(self.conf.OutputDirectory.joinpath("models"))
        os.makedirs(self.conf.OutputDirectory.joinpath("quicklook"))
        
        # Set Graphical options
        use(backend)
        self.plotformat=plotformat
        
        return None
    
    def WriteResults(self, datasets):
        """
        Write Results according to requested analysis.
        
        Parameters
        ----------
        datasets : gammapy.datasets.Datasets()
            Simulated datasets.
        """
        # Set Simuulated Datasets
        self.datasets = datasets
        
        # Add time start and stop arrays
        time_start = [d.gti.time_start.value[0] - d.gti.time_ref.value for d in self.datasets] * u.d
        time_stop  = [d.gti.time_stop.value[0]  - d.gti.time_ref.value for d in self.datasets] * u.d
        
        self.time_start = time_start.to('s')
        self.time_stop  = time_stop.to( 's')
        
        if self.conf.product=="DL4":
            if self.conf.analysis=="1D":
                self.WriteDL4InfoTable(cumulative=False)
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
                            filename_models=self.conf.OutputDirectory.joinpath("models/models.yaml")
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
        #info_table = self.datasets.info_table(cumulative=cumulative)
        info_table = utils.info_table(self.datasets, cumulative=cumulative)
        
        # Add time start and stop columns
        info_table['time_start'] = self.time_start
        info_table['time_stop' ] = self.time_stop
        
        # Write light curve table
        if cumulative:
            table_name = self.conf.OutputDirectory.joinpath('quicklook/lightcurve_cumulative.ecsv')
        else:
            table_name = self.conf.OutputDirectory.joinpath('quicklook/lightcurve.ecsv')
        
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
        e_min = self.conf.AxisEnergyReco.bounds[0].value
        e_max = self.conf.AxisEnergyReco.bounds[1].value
        
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
        figure_name = self.conf.OutputDirectory.joinpath(f"quicklook/lightcurve.{self.plotformat}")
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
        table_name = self.conf.OutputDirectory.joinpath(f"quicklook/stacked_counts.ecsv")
        self.log.info(f"Write {table_name}")
        stacked_table.write(table_name)
        
        # Write Table for spectrum
        spectrum_table = Table()
        spectrum_table['e_min'] = stacked.geoms['geom'].axes['energy'].edges_min
        spectrum_table['e_max'] = stacked.geoms['geom'].axes['energy'].edges_max
        spectrum_table['e_ref'] = stacked.geoms['geom'].axes['energy'].center
        spectrum_table['counts']= np.squeeze(stacked.counts.data)
        spectrum_table['background']= np.squeeze(stacked.background.data)
        spectrum_table['excess']= np.squeeze(stacked.excess.data)
        try:
            spectrum_table['counts_off']= np.squeeze(stacked.counts_off.data)
            spectrum_table['alpha']= np.squeeze(stacked.alpha.data)
            spectrum_table['acceptance']= np.squeeze(stacked.acceptance.data)
            spectrum_table['acceptance_off']= np.squeeze(stacked.acceptance_off.data)
            spectrum_table['npred_off']=np.squeeze(stacked.npred_off().data)
        except AttributeError as e:
            self.log.warning(f"No off quantities simulated.")
        spectrum_table['npred_signal']=np.squeeze(stacked.npred_signal().data)
        spectrum_table['npred_background']=np.squeeze(stacked.npred_background().data)
        spectrum_table['npred']=np.squeeze(stacked.npred().data)
        spectrum_table['stat_array']=np.squeeze(stacked.stat_array())
        spectrum_table['sqrt_ts']=stacked._counts_statistic[stacked.mask_safe.data].sqrt_ts
        spectrum_table['counts_rate']= spectrum_table['counts'] / stacked_dict['livetime']
        spectrum_table['background_rate']= spectrum_table['background'] / stacked_dict['livetime']
        spectrum_table['excess_rate']= spectrum_table['excess'] / stacked_dict['livetime']
        table_name = self.conf.OutputDirectory.joinpath(f"quicklook/stacked_spectrum.ecsv")
        self.log.info(f"Write {table_name}")
        spectrum_table.write(table_name)
        
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

        # Graphics
        ax.set_xlabel(f"Energy / {stacked.geoms['geom'].axes['energy'].unit}", fontsize = 'large')
        ax.set_ylabel('Counts', fontsize = 'large')
        ax.set_xscale('log')
        ax.set_yscale('log')
        
        title=f"Counts energy distribution in [{time_min:2f}, {time_max:2f}] {u.s}. T0={self.conf.timeRef}."
        ax.set_title(title, fontsize = 'large')
        ax.grid()
        ax.legend(bbox_to_anchor=(1.2, 1.0), loc="upper right")
        
        # Save Plot
        figure_name = self.conf.OutputDirectory.joinpath(f"quicklook/stacked_spectrum.{self.plotformat}")
        self.log.info(f"Write {figure_name}")
        fig.savefig(figure_name)
        return None
    
    
    def PlotCTAIRFs(self, irfs):
        """
        Plot and save the CTA IRFs.
        
        Parameters
        ----------
        irfs : dict
            Dictionary with CTA IRFs.
        """
        try:
            self.PlotCTAEffectiveArea(irfs['aeff'])
        except KeyError as e:
            self.log.warning(e)
            
        try:
            self.PlotCTABackground(irfs['bkg'])
        except KeyError as e:
            self.log.warning(e)
        
        try:
            self.PlotCTAEnergyDispersion(irfs['edisp'])
        except KeyError as e:
            self.log.warning(e)
            
        try:
            self.PlotCTAPSF(irfs['psf'])
        except KeyError as e:
            self.log.warning(e)
            
        return None
    
    def PlotCTAEffectiveArea(self, aeff):
        """
        Plot effective area.
        
        Parameters
        ----------
        aeff : gammapy.irf.EffectiveAreaTable2D
            Effective Area.
        """
        fig, axs = plt.subplots(1,3, figsize=(24,8), constrained_layout=True)
        aeff.plot_energy_dependence(ax=axs[0], offset=Angle([0,0.4,1,2,3]*u.deg), drawstyle='steps-mid')
        aeff.plot(ax=axs[1], add_cbar=True)
        aeff.plot_offset_dependence(ax=axs[2])
        
        axs[0].loglog()
        
        axs[0].grid()
        axs[1].grid(color='white', ls='dotted')
        axs[2].grid()
        
        figure_name = self.conf.OutputDirectory.joinpath(f"irfs/effectivearea.{self.plotformat}")
        self.log.info(f"Write {figure_name}")
        fig.savefig(figure_name)
        return None
    
    def PlotCTABackground(self, bkg):
        """
        Plot Background.
        
        Parameters
        ----------
        bkg : gammapy.irf.Background2D or Background3D
            Background.
        """
        return None
    
    def PlotCTAEnergyDispersion(self, edisp):
        """
        Plot Energy Dispersion.
        
        Parameters
        ----------
        edisp : gammapy.irf.
            Energy Dispersion.
        """
        return None
    
    def PlotCTAPSF(self, psf):
        """
        Plot PSF.
        
        Parameters
        ----------
        psf : gammapy.irf.
            PSF.
        """
        return None
    
    

    def PlotStep(self, functions, labels):
        """
        Plot functions as step histograms.
        
        Parameters
        ----------
        functions : list of (astropy.Quantity, astropy.Quantity, str)
            First quantity is x centroid, second is y value, third is label.
        labels : dict
            Values are string, keys are xlabel, ylabel, xscale, yscale, title, figurename.
        """
        
        # Make Plot
        fig, ax = plt.subplots(1, figsize = (10, 5), constrained_layout=True)
        
        for function in functions:
            ax.step(function[0], function[1], where='mid', label=function[2])
        
        # Graphics
        ax.set_xlabel(labels['xlabel'], fontsize = 'large')
        ax.set_ylabel(labels['ylabel'], fontsize = 'large')
        ax.set_xscale(labels['xscale'])
        ax.set_yscale(labels['yscale'])
        ax.set_title( labels['title' ], fontsize = 'large')
        ax.grid()
        ax.legend()
        
        # Save Plot
        figure_name = self.conf.OutputDirectory.joinpath(f"{labels['figurename']}.{self.plotformat}")
        self.log.info(f"Write {figure_name}")
        fig.savefig(figure_name)
        return None
    
    def PlotDRM(self,
                MapDRM,
                stretch='linear',
                cmap='plasma',
                extracolor='white',
                filename='DRM'
                ):
        """
        Plot the DRM read from file.
        
        Parameters
        ----------
        MapDRM : gammapy.maps.RegionNDMap
            Map containing the Detector Response Matrix and its geometry.
        stretch : str or float
            Normalization of the Colorbar.
        cmap :  str
            Colormap name of `matplotlib.colorbar.Colorbar`
        extracolor : str
            Color for secondary graphical elements.
        filename : str
            File name for the plot.
        """
        
        # Prepare Grid
        X, Y = np.meshgrid(MapDRM.geom.axes['energy_true'].center.value,
                           MapDRM.geom.axes['energy'].center.value
                           )

        # Copy Data with Masking
        DRM = np.squeeze(MapDRM.data)
        Z = np.ma.masked_where(DRM.T <= 0, DRM.T)
        
        # Define Color Normalization
        if stretch == "log":
            norm = LogNorm()
            
            # Define Countour Levels
            levs = np.linspace(np.floor(np.log10(Z.min())),
                               np.ceil( np.log10(Z.max())),
                               num = 50
                              )
            levs = np.power(10, levs)
        
        elif stretch in ['linear', 'sqrt'] or isinstance(stretch, float):
            
            if stretch=="linear":
                gamma=1.0
            elif stretch=="sqrt":
                gamma=0.5
            else:
                gamma=stretch
            norm = PowerNorm(gamma=gamma)
            
            # Define Countour Levels
            levs = np.linspace(np.floor(np.power(Z.min(),gamma)),
                               np.ceil( np.power(Z.max(),gamma)),
                               num = 50
                              )
            levs = np.power(levs, 1.0/gamma)
        else:
            raise NotImplementedError(f"stretch={stretch} not allowed. Only float or string in[\'linear\', \'sqrt\', \'log\'] are allowed.") 


        # Plot
        fig, ax = plt.subplots(1, figsize=(9,5))

        # Plot Data
        cs = ax.contourf(X, Y, Z, levs, norm = norm, cmap = cmap)
        _  = ax.contour( X, Y, Z, levs, norm = norm, colors=extracolor, alpha=0.05)
        cbar = fig.colorbar(cs, ax=ax, location='right', shrink=0.9)
        
        # Graphics
        ax.set_facecolor(cs.cmap(0))
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(color=extracolor, ls='dotted')
        cbar.set_label(f"Effective Area Redistribution / {MapDRM.unit}", fontsize = 'large', rotation=90)
        ax.set_xlabel(f"True Energy [{MapDRM.geom.axes['energy_true'].unit}]", fontsize = 'large')
        ax.set_ylabel(f"Reco Energy [{MapDRM.geom.axes['energy'].unit}]", fontsize = 'large')
        ax.set_title(f"Detector Response Matrix {self.conf.instrument} {self.conf.detector}",
                     fontsize = 'large'
                     )

        # Save Plot
        figure_name = self.conf.OutputDirectory.joinpath(f"irfs/{filename}.{self.plotformat}")
        self.log.info(f"Write {figure_name}")
        fig.savefig(figure_name, facecolor = 'white')
        
        return None
    
    
    def PlotTemporalModel(self, temporal_model):
        """
        Plot and save the temporal model.
        
        Parameters
        ----------
        temporal_model : gammapy.modelling.models.TemporalModel
            Temporal Model to plot.
        """
        # Plot
        fig, ax = plt.subplots(1, figsize=(9,5))
        
        PlotRange= (self.conf.timeRef+self.conf.timeStart, self.conf.timeRef+self.conf.timeStop)
        temporal_model.plot(PlotRange, ax=ax)
        
        # Save Plot
        figure_name = self.conf.OutputDirectory.joinpath(f"models/temporal.{self.plotformat}")
        self.log.info(f"Write {figure_name}")
        fig.savefig(figure_name, facecolor = 'white')
        
        return None
