import pandas as pd
import h5py as h5
from cosmic.evolve import Evolve
from cosmic._version import __version__
from cosmic.plotting import plot_binary_evol
import matplotlib.pyplot as plt
import warnings


__all__ = ['COSMICOutput', 'save_initC', 'load_initC']


class COSMICOutput:
    def __init__(self, bpp=None, bcm=None, initC=None, kick_info=None, file=None, label=None):
        """Container for COSMIC output data components.

        Can be initialized either from data components directly or by loading from an HDF5 file.

        Parameters
        ----------
        bpp : `pandas.DataFrame`, optional
            Important evolution timestep table, by default None
        bcm : `pandas.DataFrame`, optional
            User-defined timestep table, by default None
        initC : `pandas.DataFrame`, optional
            Initial conditions table, by default None
        kick_info : `pandas.DataFrame`, optional
            Natal kick information table, by default None
        file : `str`, optional
            Filename/path to HDF5 file to load data from, by default None
        label : `str`, optional
            Optional label for the output instance, by default None

        Raises
        ------
        ValueError
            If neither file nor all data components are provided.
        """
        # require that either file is given or all data components are given
        if file is None and (bpp is None or bcm is None or initC is None or kick_info is None):
            raise ValueError("Either file or all data components (bpp, bcm, initC, kick_info) must be provided.")
        if file is not None:
            self.bpp = pd.read_hdf(file, key='bpp')
            self.bcm = pd.read_hdf(file, key='bcm')
            self.initC = load_initC(file, key='initC', settings_key='initC_settings')
            self.kick_info = pd.read_hdf(file, key='kick_info')
            with h5.File(file, 'r') as f:
                file_version = f.attrs.get('COSMIC_version', 'unknown')
                label = f.attrs.get('label', '')
            self.label = label if label != '' else None
            if file_version != __version__:
                warnings.warn(f"You have loaded COSMICOutput from a file that was run using COSMIC version {file_version}, "
                              f"but the current version is {__version__}. "
                              "There may be compatibility issues, or differences in output when rerunning, be sure to check the changelog.", UserWarning)
        else:
            self.bpp = bpp
            self.bcm = bcm
            self.initC = initC
            self.kick_info = kick_info
            self.label = label if label is not None else None

    def __len__(self):
        return len(self.initC)

    def __repr__(self):
        return f'<COSMICOutput{" - " + self.label if self.label is not None else ""}: {len(self)} {"binaries" if len(self) != 1 else "binary"}>'
        
    @property
    def final_bpp(self):
        """Get the final timestep for each binary from the bpp table.

        Returns
        -------
        final_bpp : `pandas.DataFrame`
            DataFrame containing only the final timestep for each binary.
        """
        return self.bpp.drop_duplicates(subset='bin_num', keep='last')

    def save(self, output_file):
        """Save all data components to an HDF5 file

        Parameters
        ----------
        output_file : `str`
            Filename/path to the HDF5 file
        """
        self.bpp.to_hdf(output_file, key='bpp')
        self.bcm.to_hdf(output_file, key='bcm')
        save_initC(output_file, self.initC, key='initC', settings_key='initC_settings')
        self.kick_info.to_hdf(output_file, key='kick_info')
        with h5.File(output_file, 'a') as f:
            f.attrs['COSMIC_version'] = __version__
            f.attrs['label'] = self.label if self.label is not None else ''

    def rerun_with_settings(self, new_settings, reset_kicks=False, inplace=False):
        """Rerun the simulation with new settings.

        Parameters
        ----------
        new_settings : `dict`
            Dictionary of new settings to apply. Any setting not included will retain its original value.
        reset_kicks : `bool`, optional
            If True, reset natal kicks to be randomly sampled again.
            If False, retain original kicks. By default False.
            (You may want to reset the kicks if changing settings that affect remnant masses or
            kick distribution.)
        inplace : `bool`, optional
            If True, update the current instance. If False, return a new instance. By default False.
        
        Returns
        -------
        new_output : `COSMICOutput`
            New COSMICOutput instance with updated simulation results (only if inplace is False).
        """
        # merge new settings with existing initC
        updated_initC = self.initC.copy()
        for key, value in new_settings.items():
            if key in updated_initC.columns:
                updated_initC[key] = value
            else:
                raise KeyError(f"Setting '{key}' not found in initC columns.")
            
        # reset kicks if requested
        if reset_kicks:
            kick_cols = ["natal_kick_1", "natal_kick_2", "phi_1", "phi_2", "theta_1", "theta_2",
                         "mean_anomaly_1", "mean_anomaly_2"]
            for col in kick_cols:
                updated_initC[col] = -100.0
        elif 'kickflag' in new_settings or 'remnantflag' in new_settings:
            warnings.warn(
                "You have changed 'kickflag' or 'remnantflag' without resetting kicks. "
                "This may lead to inconsistent results if the kick distribution or remnant masses have changed. "
                "Consider setting reset_kicks=True.", UserWarning
            )

        # re-run the simulation
        new_bpp, new_bcm, new_initC, new_kick_info = Evolve.evolve(initialbinarytable=updated_initC)
        
        if inplace:
            self.bpp = new_bpp
            self.bcm = new_bcm
            self.initC = new_initC
            self.kick_info = new_kick_info
        else:
            return COSMICOutput(bpp=new_bpp, bcm=new_bcm, initC=new_initC, kick_info=new_kick_info)
        

    def plot_detailed_evolution(self, bin_num, show=True, **kwargs):
        """Plot detailed evolution for a specific binary.

        Parameters
        ----------
        bin_num : `int`
            Index of the binary to plot.
        **kwargs : 
            Additional keyword arguments passed to the plotting function (plotting.plot_binary_evol).
        """
        # check the bin_num is in the bcm
        if bin_num not in self.bcm['bin_num'].values:
            raise ValueError(f"bin_num {bin_num} not found in bcm table.")

        # warn if bcm has only two entries for this binary
        bcm_subset = self.bcm[self.bcm['bin_num'] == bin_num]
        if len(bcm_subset) <= 2:
            warnings.warn(
                f"bcm table for bin_num {bin_num} has only {len(bcm_subset)} entries. Detailed evolution "
                "plot may be uninformative. You should set dtp, or timestep_conditions, to increase the "
                "number of timesteps in the bcm table.", UserWarning
            )

        fig = plot_binary_evol(self.bcm.loc[bin_num], **kwargs)
        if show:
            plt.show()
        return fig


    def plot_distribution(self, x_col, y_col=None, c_col=None, when='final',
                          fig=None, ax=None, show=True,
                          xlabel='auto', ylabel='auto', clabel='auto', **kwargs):
        """Plot distribution of binaries in specified columns.

        Plots can be histograms (if only x_col is given) or scatter plots (if both x_col and y_col are given).
        Optionally, colour coding can be applied using c_col.

        Parameters
        ----------
        x_col : `str`
            Column name for x-axis.
        y_col : `str`, optional
            Column name for y-axis. If None, a histogram will be plotted. By default None.
        c_col : `str`, optional
            Column name for colour coding. By default None.
        when : `str`, optional
            When to take the values from: 'initial' or 'final'. By default 'final'.
        fig : `matplotlib.figure.Figure`, optional
            Figure to plot on. If None, a new figure is created. By default None.
        ax : `matplotlib.axes.Axes`, optional
            Axes to plot on. If None, new axes are created. By default None.
        show : `bool`, optional
            If True, display the plot immediately. By default True.
        xlabel : `str`, optional
            Label for x-axis. If 'auto', uses the column name. By default 'auto'.
        ylabel : `str`, optional
            Label for y-axis. If 'auto', uses the column name or 'Count' for histogram. By default 'auto'.
        clabel : `str`, optional
            Label for colorbar. If 'auto', uses the column name. By default 'auto
        **kwargs :
            Additional keyword arguments passed to the plotting function.

        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The figure containing the plot.
        ax : `matplotlib.axes.Axes`
            The axes containing the plot.
        """
        if fig is None or ax is None:
            fig, ax = plt.subplots()

        if when == 'initial':
            data = self.initC
        elif when == 'final':
            data = self.bpp.drop_duplicates(subset='bin_num', keep='last')
        else:
            raise ValueError("Parameter 'when' must be either 'initial' or 'final'.")
        
        if xlabel == 'auto':
            xlabel = x_col
        if ylabel == 'auto':
            ylabel = y_col if y_col is not None else 'Count'
        if clabel == 'auto' and c_col is not None:
            clabel = c_col
        
        if y_col is None:
            # histogram
            ax.hist(data[x_col], bins=kwargs.get('bins', "fd"),
                    color=kwargs.get('color', "tab:blue"), **kwargs)
            ax.set(
                xlabel=xlabel,
                ylabel=ylabel,
            )
        else:
            # scatter plot
            sc = ax.scatter(data[x_col], data[y_col],
                            c=data[c_col] if c_col is not None else kwargs.get('color', "tab:blue"),
                            **kwargs)
            ax.set(
                xlabel=xlabel,
                ylabel=ylabel,
            )
            if c_col is not None:
                cbar = fig.colorbar(sc, ax=ax)
                cbar.set_label(clabel)

        if show:
            plt.show()
        return fig, ax


def save_initC(filename, initC, key="initC", settings_key="initC_settings", force_save_all=False):
    """Save an initC table to an HDF5 file.

    Any column where every binary has the same value (setting) is saved separately with only a single copy
    to save space.

    This will take slightly longer (a few seconds instead of 1 second) to run but will save you around
    a kilobyte per binary, which adds up!

    Parameters
    ----------
    filename : `str`
        Filename/path to the HDF5 file
    initC : `pandas.DataFrame`
        Initial conditions table
    key : `str`, optional
        Dataset key to use for main table, by default "initC"
    settings_key : `str`, optional
        Dataset key to use for settings table, by default "initC_settings"
    force_save_all : `bool`, optional
        If true, force all settings columns to be saved in the main table, by default False
    """

    # for each column, check if all values are the same
    uniques = initC.nunique(axis=0)
    compress_cols = [col for col in initC.columns if uniques[col] == 1]

    if len(compress_cols) == 0 or force_save_all:
        # nothing to compress, just save the whole table
        initC.to_hdf(filename, key=key)
    else:
        # save the main table without the compressed columns
        initC.drop(columns=compress_cols).to_hdf(filename, key=key)

        # save the compressed columns separately
        settings_df = pd.DataFrame([{col: initC[col].iloc[0] for col in compress_cols}])
        settings_df.to_hdf(filename, key=settings_key)


def load_initC(filename, key="initC", settings_key="initC_settings"):
    """Load an initC table from an HDF5 file.

    If settings were saved separately, they are merged back into the main table.

    Parameters
    ----------
    filename : `str`
        Filename/path to the HDF5 file
    key : `str`, optional
        Dataset key to use for main table, by default "initC"
    settings_key : `str`, optional
        Dataset key to use for settings table, by default "initC_settings"

    Returns
    -------
    initC : `pandas.DataFrame`
        Initial conditions table
    """

    with h5.File(filename, 'r') as f:
        has_settings = settings_key in f.keys()

    initC = pd.read_hdf(filename, key=key)

    if has_settings:
        settings_df = pd.read_hdf(filename, key=settings_key)
        initC.loc[:, settings_df.columns] = settings_df.values[0]

    return initC


