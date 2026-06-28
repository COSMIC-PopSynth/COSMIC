"""Vectorized parameter space definition.

A `ParameterSpace` holds an ordered list of `Parameter` instances and provides
vectorized operations on (N, D) sample arrays including sampling, scale
transforms, prior evaluation, and bounds checking.
"""
import numpy as np
from dataclasses import dataclass

from .distributions import Distribution, get_distribution


@dataclass
class Parameter:
    """A single dimension in the parameter space.

    Parameters
    ----------
    name : `str`
        Name of the parameter (used for column ordering).
    min_value : `float`
        Lower bound, always in physical space (e.g. solar masses for
        ``'kroupa'``, days for ``'sana'``, km/s for ``'log_normal'``).  The
        parameter's distribution maps this into sampling space via its
        transform.
    max_value : `float`
        Upper bound, in physical space (same convention as ``min_value``).
    dist : `str` or `~cosmic.sample.stroopwafel.distributions.Distribution`, optional
        The prior distribution, given either as a name registered in
        :data:`~cosmic.sample.stroopwafel.distributions.DISTRIBUTIONS`
        (e.g. ``'kroupa'``, ``'sana'``, ``'flat_in_log'``) or as a
        :class:`~cosmic.sample.stroopwafel.distributions.Distribution`
        instance for custom priors.  By default ``'uniform'``.

    Attributes
    ----------
    distribution : `Distribution`
        The resolved distribution instance.
    lo, hi : `float`
        The bounds in sampling space (``min_value``/``max_value`` passed
        through ``distribution.transform``).
    """
    name: str
    min_value: float
    max_value: float
    dist: "str | Distribution" = 'uniform'

    def __post_init__(self):
        self.distribution = get_distribution(self.dist)
        self.lo, self.hi = self.distribution.transform.bounds(self.min_value, self.max_value)


class ParameterSpace:
    """An ordered collection of Parameters with vectorized operations.

    All methods operate on (N, D) numpy arrays where columns are ordered
    alphabetically by parameter name.

    Parameters
    ----------
    params : `list` of `Parameter`
        List of parameter definitions. They will be sorted by name
        internally.
    """

    def __init__(self, params):
        # Sort by name for deterministic column ordering (same as old code)
        self.params = sorted(params, key=lambda p: p.name)
        self.names = [p.name for p in self.params]
        self._name_to_idx = {p.name: i for i, p in enumerate(self.params)}
        self.ndim = len(self.params)

    def idx(self, name):
        """Get the column index for a parameter by name.

        Parameters
        ----------
        name : `str`
            Parameter name.

        Returns
        -------
        `int`
            Column index in the (N, D) sample arrays.
        """
        return self._name_to_idx[name]

    def sample(self, n, rng=None):
        """Draw n samples from the prior distribution.

        Parameters
        ----------
        n : `int`
            Number of samples to draw.
        rng : `numpy.random.Generator`, optional
            Random number generator, by default None

        Returns
        -------
        samples : `numpy.ndarray`
            (N, D) array of samples in sampling space.
        mask : `numpy.ndarray`
            (N,) boolean array indicating which samples are in bounds.
        """
        rng = rng or np.random.default_rng()
        samples = np.empty((n, self.ndim))
        mask = np.ones(n, dtype=bool)

        for i, p in enumerate(self.params):
            col = p.distribution.sample(n, p.lo, p.hi, rng=rng)
            samples[:, i] = col
            mask &= (col >= p.lo) & (col <= p.hi)

        return samples, mask

    def in_bounds(self, samples):
        """Check which rows of an (N, D) array are within parameter bounds.

        Parameters
        ----------
        samples : `numpy.ndarray`
            (N, D) array of samples in sampling space.

        Returns
        -------
        `numpy.ndarray`
            (N,) boolean mask where True means the sample is in bounds.
        """
        mask = np.ones(len(samples), dtype=bool)
        for i, p in enumerate(self.params):
            mask &= (samples[:, i] >= p.lo) & (samples[:, i] <= p.hi)
        return mask

    def to_physical(self, samples):
        """Convert an (N, D) array from sampling space to physical space.

        Parameters
        ----------
        samples : `numpy.ndarray`
            (N, D) array in sampling space.

        Returns
        -------
        `numpy.ndarray`
            (N, D) array in physical space.
        """
        result = samples.copy()
        for i, p in enumerate(self.params):
            result[:, i] = p.distribution.transform.to_physical(samples[:, i])
        return result

    def to_sampling(self, samples):
        """Convert an (N, D) array from physical space to sampling space.

        Parameters
        ----------
        samples : `numpy.ndarray`
            (N, D) array in physical space.

        Returns
        -------
        `numpy.ndarray`
            (N, D) array in sampling space.
        """
        result = samples.copy()
        for i, p in enumerate(self.params):
            result[:, i] = p.distribution.transform.to_sampling(samples[:, i])
        return result

    def compute_prior(self, samples):
        """Compute the prior probability for each row.

        The joint prior is the product of the marginal priors across all
        dimensions.

        Parameters
        ----------
        samples : `numpy.ndarray`
            (N, D) array in sampling space.

        Returns
        -------
        `numpy.ndarray`
            (N,) array of prior probabilities.
        """
        log_prior = np.zeros(len(samples))
        for i, p in enumerate(self.params):
            col_prior = p.distribution.pdf(samples[:, i], p.lo, p.hi)
            # Clamp to avoid log(0)
            log_prior += np.log(np.maximum(col_prior, 1e-300))
        return np.exp(log_prior)

    def compute_sigma(self, hit_samples, average_density_one_dim):
        """Compute per-hit, per-dimension Gaussian widths (sigma).

        Parameters
        ----------
        hit_samples : `numpy.ndarray`
            (K, D) array of hit locations in sampling space.
        average_density_one_dim : `float`
            Characteristic inter-sample spacing, typically
            ``1 / num_explored ** (1 / D)``.

        Returns
        -------
        `numpy.ndarray`
            (K, D) array of sigma values.
        """
        K = len(hit_samples)
        sigmas = np.empty((K, self.ndim))

        for i, p in enumerate(self.params):
            try:
                sigmas[:, i] = p.distribution.sigma(
                    hit_samples[:, i], p.lo, p.hi, average_density_one_dim
                )
            except ValueError as e:
                raise ValueError(f"Parameter '{p.name}': {e}") from e
        return sigmas
