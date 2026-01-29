import json
import pandas as pd
import h5py as h5
from cosmic.evolve import Evolve
from cosmic._version import __version__
from cosmic.plotting import plot_binary_evol
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np
import warnings


__all__ = ['COSMICOutput', 'COSMICPopOutput', 'save_initC', 'load_initC']


kstar_translator = [
    {'long': 'Main Sequence (Low mass)', 'short': 'MS < 0.7', 'colour': (0.996078, 0.843476, 0.469158, 1.0)},
    {'long': 'Main Sequence', 'short': 'MS', 'colour': (0.996078, 0.843476, 0.469158, 1.0)},
    {'long': 'Hertzsprung Gap', 'short': 'HG', 'colour': (0.939608, 0.471373, 0.094902, 1.0)},
    {'long': 'First Giant Branch', 'short': 'FGB', 'colour': (0.716186, 0.833203, 0.916155, 1.0)},
    {'long': 'Core Helium Burning', 'short': 'CHeB', 'colour': (0.29098, 0.59451, 0.78902, 1.0)},
    {'long': 'Early AGB', 'short': 'EAGB', 'colour': (0.294902, 0.690196, 0.384314, 1.0)},
    {'long': 'Thermally Pulsing AGB', 'short': 'TPAGB',
     'colour': (0.723122, 0.889612, 0.697178, 1.0)},
    {'long': 'Helium Main Sequence', 'short': 'HeMS', 'colour': (0.254627, 0.013882, 0.615419, 1.0)},
    {'long': 'Helium Hertsprung Gap', 'short': 'HeHG', 'colour': (0.562738, 0.051545, 0.641509, 1.0)},
    {'long': 'Helium Giant Branch', 'short': 'HeGB', 'colour': (0.798216, 0.280197, 0.469538, 1.0)},
    {'long': 'Helium White Dwarf', 'short': 'HeWD', 'colour': (0.368166, 0.232828, 0.148275, 1.0)},
    {'long': 'Carbon/Oxygen White Dwarf', 'short': 'COWD', 'colour': (0.620069, 0.392132, 0.249725, 1.0)},
    {'long': 'Oxygen/Neon White Dwarf', 'short': 'ONeWD', 'colour': (0.867128, 0.548372, 0.349225, 1.0)},
    {'long': 'Neutron Star', 'short': 'NS', 'colour': (0.501961, 0.501961, 0.501961, 1.0)},
    {'long': 'Black Hole', 'short': 'BH', 'colour': (0.0, 0.0, 0.0, 1.0)},
    {'long': 'Massless Remnant', 'short': 'MR', 'colour': "white"},
    {'long': 'Chemically Homogeneous', 'short': 'CHE', 'colour': (0.647059, 0.164706, 0.164706, 1.0)}
]


class COSMICOutput:
    def __init__(self, bpp=None, bcm=None, initC=None, kick_info=None, file=None, label=None,
                 file_key_suffix=''):
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
        file_key_suffix : `str`, optional
            Suffix to append to dataset keys when loading from file, by default ''. E.g. if set to '_singles',
            datasets 'bpp_singles', 'bcm_singles', etc. will be loaded as bpp, bcm, etc.

        Raises
        ------
        ValueError
            If neither file nor all data components are provided.
        """
        # require that either file is given or all data components are given
        if file is None and (bpp is None or bcm is None or initC is None or kick_info is None):
            raise ValueError("Either file or all data components (bpp, bcm, initC, kick_info) must be provided.")
        if file is not None:
            self.bpp = pd.read_hdf(file, key=f'bpp{file_key_suffix}')
            self.bcm = pd.read_hdf(file, key=f'bcm{file_key_suffix}')
            self.initC = load_initC(file, key=f'initC{file_key_suffix}',
                                    settings_key=f'initC_{file_key_suffix}_settings')
            self.kick_info = pd.read_hdf(file, key=f'kick_info{file_key_suffix}')
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
    
    def __getitem__(self, key):
        """Subselect binaries by bin_num across all data components.
        Keys can be integers or lists/arrays of integers or slices.
        If the key is an array of bools, mask initC to get the corresponding bin_nums."""
        # convert key to list of bin_nums, regardless of input type
        if isinstance(key, int):
            key = [key]
        elif isinstance(key, slice):
            key = self.initC['bin_num'].iloc[key].tolist()
        elif isinstance(key, (pd.Series, list, np.ndarray)) and len(key) == len(self.initC) and isinstance(key[0], (bool, np.bool_)):
            if not key.any():
                raise IndexError("Boolean mask resulted in zero selected binaries.")
            key = self.initC['bin_num'][key].tolist()
        # otherwise, reject invalid types
        elif not isinstance(key, (list, np.ndarray, pd.Series)):
            raise TypeError("Key must be an int, slice, list/array of ints, or boolean mask.")

        bpp_subset = self.bpp[self.bpp['bin_num'].isin(key)]
        bcm_subset = self.bcm[self.bcm['bin_num'].isin(key)]
        initC_subset = self.initC[self.initC['bin_num'].isin(key)]
        kick_info_subset = self.kick_info[self.kick_info['bin_num'].isin(key)]
        return COSMICOutput(bpp=bpp_subset, bcm=bcm_subset, initC=initC_subset, kick_info=kick_info_subset, label=self.label)
        
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

        if "ktype_kwargs" not in kwargs:
            kwargs["ktype_kwargs"] = {'k_type_colors': [kstar_translator[k]["colour"] for k in range(len(kstar_translator))]}
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
            c = data[c_col] if c_col is not None else kwargs.get('color', None)
            if c_col == 'kstar_1' or c_col == 'kstar_2':
                c = data[c_col].map(lambda k: kstar_translator[k]['colour'])
            sc = ax.scatter(data[x_col], data[y_col],
                            c=c,
                            **kwargs)
            ax.set(
                xlabel=xlabel,
                ylabel=ylabel,
            )
            if c_col is not None and c_col not in ['kstar_1', 'kstar_2']:
                cbar = fig.colorbar(sc, ax=ax)
                cbar.set_label(clabel)
            elif c_col is not None:
                # extract colours and labels
                colours = [entry["colour"] for entry in kstar_translator[1:-2]]
                labels = [entry["short"] for entry in kstar_translator[1:-2]]

                # create colormap
                cmap = ListedColormap(colours)
                bounds = np.arange(len(colours) + 1)
                norm = BoundaryNorm(bounds, cmap.N)

                cb = plt.colorbar(
                    mappable=plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                    ticks=np.arange(len(colours)) + 0.5,
                    boundaries=bounds,
                    ax=ax
                )
                cb.ax.set_yticklabels(labels)
                cb.set_label(clabel)

        if show:
            plt.show()
        return fig, ax
    

class COSMICPopOutput():
    def __init__(self, file, label=None):
        # read in convergence tables and totals
        keys = ['conv', 'idx', 'match', 'mass_binaries', 'mass_singles',
                'n_binaries', 'n_singles', 'mass_stars', 'n_stars']
        for key in keys:
            setattr(self, key, pd.read_hdf(file, key=key))

        # load config back from JSON storage
        with h5.File(file, 'r') as f:
            self.config = json.loads(f['config'][()])

        # create a COSMICOutput for the binaries, and optionally for the singles
        self.output = COSMICOutput(file=file, label=label + ' [binaries]' if label is not None else None)
        singles = "keep_singles" in self.config["sampling"] and self.config["sampling"]["keep_singles"]
        self.singles_output = COSMICOutput(
            file=file, label=label + ' [singles]' if label is not None else None,
            file_key_suffix='_singles'
        ) if singles else None
        self.label = label

    def __repr__(self):
        r = f'<COSMICPopOutput{" - " + self.label if self.label is not None else ""}: {len(self.output)} binaries>'
        if self.singles_output is not None:
            r = r[:-1] + f', {len(self.singles_output)} singles>'
        return r
    
    def __len__(self):
        return len(self.conv)
    
    def to_combined_output(self):
        """Combine binaries and singles into a single COSMICOutput instance.

        Returns
        -------
        combined_output : `COSMICOutput`
            COSMICOutput instance containing both binaries and singles.
        
        Raises
        ------
        ValueError
            If singles output is not available.
        """
        if self.singles_output is None:
            raise ValueError("Singles output is not available in this COSMICPopOutput instance.")

        bpp = pd.concat([self.output.bpp, self.singles_output.bpp], ignore_index=True)
        bcm = pd.concat([self.output.bcm, self.singles_output.bcm], ignore_index=True)
        initC = pd.concat([self.output.initC, self.singles_output.initC], ignore_index=True)
        kick_info = pd.concat([self.output.kick_info, self.singles_output.kick_info], ignore_index=True)

        return COSMICOutput(
            bpp=bpp,
            bcm=bcm,
            initC=initC,
            kick_info=kick_info,
            label=self.label + ' [combined]' if self.label is not None else None
        )


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


