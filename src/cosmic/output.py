import pandas as pd
import h5py as h5
from cosmic.evolve import Evolve


__all__ = ['COSMICOutput', 'save_initC', 'load_initC']


class COSMICOutput:
    def __init__(self, bpp=None, bcm=None, initC=None, kick_info=None, file=None):
        # require that either file is given or all data components are given
        if file is None and (bpp is None or bcm is None or initC is None or kick_info is None):
            raise ValueError("Either file or all data components (bpp, bcm, initC, kick_info) must be provided.")
        if file is not None:
            self.bpp = pd.read_hdf(file, key='bpp')
            self.bcm = pd.read_hdf(file, key='bcm')
            self.initC = load_initC(file, key='initC', settings_key='initC_settings')
            self.kick_info = pd.read_hdf(file, key='kick_info')
        else:
            self.bpp = bpp
            self.bcm = bcm
            self.initC = initC
            self.kick_info = kick_info

    def __len__(self):
        return len(self.initC)

    def __repr__(self):
        return f'<COSMICOutput: {len(self)} {"binaries" if len(self) != 1 else "binary"}>'
        
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

    def rerun_with_settings(self, new_settings, inplace=False):
        """Rerun the simulation with new settings.

        Parameters
        ----------
        new_settings : `dict`
            Dictionary of new settings to apply. Any setting not included will retain its original value.
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

        # Rerun the simulation
        new_bpp, new_bcm, new_initC, new_kick_info = Evolve.evolve(initialbinarytable=updated_initC)
        
        if inplace:
            self.bpp = new_bpp
            self.bcm = new_bcm
            self.initC = new_initC
            self.kick_info = new_kick_info
        else:
            return COSMICOutput(bpp=new_bpp, bcm=new_bcm, initC=new_initC, kick_info=new_kick_info)


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


