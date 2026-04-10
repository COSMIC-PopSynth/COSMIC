"""HDF5-based I/O for STROOPWAFEL state and results.

Replaces the old CSV-based print_samples/read_samples with array-native storage.
"""
import pandas as pd
import h5py


def save_result(path, result):
    """Save a `STROOPWAFELResult` to an HDF5 file.

    Parameters
    ----------
    path : `str`
        Output file path (created or overwritten).
    result : `STROOPWAFELResult`
        Result container to serialise.
    """
    with h5py.File(path, 'w') as f:
        f.create_dataset('samples', data=result.samples)
        f.create_dataset('weights', data=result.weights)
        f.create_dataset('is_hit', data=result.is_hit)
        f.create_dataset('generation', data=result.generation)
        f.create_dataset('gaussian_idx', data=result.gaussian_idx)
        f.attrs['param_names'] = result.param_names
        f.attrs['num_explored'] = result.num_explored
        f.attrs['num_hits'] = result.num_hits
        f.attrs['fraction_explored'] = result.fraction_explored


def save_cosmic_output(path, bpp_frames, initC_frames, kick_info_frames):
    """Append COSMIC output DataFrames to an HDF5 file.

    Parameters
    ----------
    path : `str`
        Output file path (appended to if it already exists).
    bpp_frames : `list` of `pandas.DataFrame`
        Binary population parameter tables from each batch.
    initC_frames : `list` of `pandas.DataFrame`
        Initial conditions tables from each batch.
    kick_info_frames : `list` of `pandas.DataFrame`
        Natal kick information tables from each batch.
    """
    if bpp_frames:
        full_bpp = pd.concat(bpp_frames, ignore_index=True)
        full_bpp.to_hdf(path, key='bpp', mode='a', format='table')
    if initC_frames:
        full_initC = pd.concat(initC_frames, ignore_index=True)
        full_initC.to_hdf(path, key='initC', mode='a', format='table')
    if kick_info_frames:
        full_kicks = pd.concat(kick_info_frames, ignore_index=True)
        full_kicks.to_hdf(path, key='kick_info', mode='a', format='table')


def save_mixture(path, mixture, generation):
    """Save `GaussianMixture` state to an HDF5 file.

    Parameters
    ----------
    path : `str`
        Output file path (appended to if it already exists).
    mixture : `GaussianMixture`
        Mixture model to serialise.
    generation : `int`
        Refinement generation number, used as the HDF5 group key.
    """
    with h5py.File(path, 'a') as f:
        grp_name = f'mixture/gen_{generation}'
        if grp_name in f:
            del f[grp_name]
        grp = f.create_group(grp_name)
        grp.create_dataset('means', data=mixture.means)
        grp.create_dataset('covariances', data=mixture.covariances)
        grp.create_dataset('alphas', data=mixture.alphas)
        grp.attrs['rejection_rate'] = mixture.rejection_rate


def load_mixture(path, generation):
    """Load `GaussianMixture` state from an HDF5 file.

    Parameters
    ----------
    path : `str`
        File path to read from.
    generation : `int`
        Refinement generation number to load.

    Returns
    -------
    `GaussianMixture` or None
        The stored mixture model, or None if the requested generation
        is not found or the file does not exist.
    """
    from .mixture_model import GaussianMixture

    try:
        with h5py.File(path, 'r') as f:
            grp = f[f'mixture/gen_{generation}']
            return GaussianMixture(
                means=grp['means'][:],
                covariances=grp['covariances'][:],
                alphas=grp['alphas'][:],
                rejection_rate=grp.attrs['rejection_rate'],
            )
    except (KeyError, FileNotFoundError):
        return None
