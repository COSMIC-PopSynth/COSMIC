"""Vectorized parameter space definition.

A `ParameterSpace` holds an ordered list of `Parameter` instances and provides
vectorized operations on (N, D) sample arrays including sampling, scale
transforms, prior evaluation, and bounds checking.
"""
import numpy as np
from dataclasses import dataclass

from .samplers import SAMPLERS
from .priors import PRIORS
from .transforms import to_sampling_space, to_physical_space, transform_bounds
from .constants import ALPHA_IMF, SANA_G, SANA_ECC


@dataclass
class Parameter:
    """A single dimension in the parameter space.

    Parameters
    ----------
    name : `str`
        Name of the parameter (used for column ordering).
    min_value : `float`
        Lower bound.  For most samplers this is in physical space (e.g.
        solar masses for ``'kroupa'``, eccentricity for ``'sana_ecc'``).

        **Exception:** ``'sana'`` (orbital period) operates in
        log10(period / days) internally, so ``min_value`` and
        ``max_value`` must be passed in log10 space.  For example, to
        span periods from ~1.4 d to ~316 000 d use
        ``min_value=0.15, max_value=5.5`` (i.e. log10 of those values).
    max_value : `float`
        Upper bound (same unit convention as ``min_value``).
    sampler : `str`, optional
        Name of the sampling distribution, by default ``'uniform'``
    prior : `str`, optional
        Name of the prior distribution, by default ``'uniform'``
    """
    name: str
    min_value: float
    max_value: float
    sampler: str = 'uniform'
    prior: str = 'uniform'

    def __post_init__(self):
        self.lo, self.hi = transform_bounds(self.min_value, self.max_value, self.sampler)


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
            sampler_fn = SAMPLERS[p.sampler]
            col = sampler_fn(n, p.lo, p.hi, rng=rng)
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
            result[:, i] = to_physical_space(samples[:, i], p.sampler)
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
            result[:, i] = to_sampling_space(samples[:, i], p.sampler)
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
            prior_fn = PRIORS[p.prior]
            col_prior = prior_fn(samples[:, i], p.lo, p.hi)
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
            col = hit_samples[:, i]
            if p.sampler in ('kroupa', 'sana', 'sana_ecc'):
                sigmas[:, i] = self._sigma_power_law(p, col, average_density_one_dim)
            else:
                # uniform, flat_in_log, etc: sigma = avg_density / prior
                prior_fn = PRIORS[p.prior]
                prior_vals = prior_fn(col, p.lo, p.hi)
                sigmas[:, i] = average_density_one_dim / prior_vals
        return sigmas

    def _sigma_power_law(self, param, values, avg_density):
        """Compute sigma for power-law distributions.

        Maps each hit value to its CDF position, steps by ``avg_density``
        in CDF space, maps back to parameter space, and returns the
        maximum of the two resulting distances.

        The CDF of a power-law p(x) ∝ x^α on [lo, hi] is

            F(x) = (x^a − lo^a) / (hi^a − lo^a),   a = α + 1

        so the inverse is F⁻¹(u) = (u·(hi^a − lo^a) + lo^a)^{1/a}.
        All intermediate quantities stay in [lo^a, hi^a], avoiding the
        catastrophic cancellation that occurred in the previous
        ``inv_lo``-based formulation when ``lo`` was very small.

        Parameters
        ----------
        param : `Parameter`
            The parameter definition for this dimension.
        values : `numpy.ndarray`
            (K,) array of hit values in sampling space.
        avg_density : `float`
            Characteristic inter-sample spacing.

        Returns
        -------
        `numpy.ndarray`
            (K,) array of sigma values.

        Raises
        ------
        `ValueError`
            If ``param.lo <= 0`` (the CDF formula requires a positive
            lower bound) or if the sampler is not a recognised power law.
        """
        if param.sampler == 'kroupa':
            alpha = ALPHA_IMF
        elif param.sampler == 'sana':
            alpha = SANA_G
        elif param.sampler == 'sana_ecc':
            alpha = SANA_ECC
        else:
            raise ValueError(f"Unknown power-law sampler: {param.sampler}")

        if param.lo <= 0:
            raise ValueError(
                f"Parameter '{param.name}' has lo={param.lo} <= 0 but uses "
                f"the '{param.sampler}' power-law sampler.  The CDF formula "
                f"requires lo > 0.  For 'sana' (period), pass the log10(P) "
                f"lower bound directly, e.g. min_value=0.15 rather than "
                f"min_value=np.log10(0.15)."
            )

        a = alpha + 1   # 0.55 for sana_ecc, 0.45 for sana, -1.3 for kroupa
        lo_a = float(param.lo) ** a
        hi_a = float(param.hi) ** a
        range_a = hi_a - lo_a  # always finite; no huge intermediates

        # --- Forward CDF: F(x) = (x^a - lo_a) / range_a ----------------
        u = np.clip((np.power(values, a) - lo_a) / range_a, 0.0, 1.0)

        # --- Step in CDF space ------------------------------------------
        u_right = np.clip(u + avg_density, 0.0, 1.0)
        u_left  = np.clip(u - avg_density, 0.0, 1.0)

        # --- Inverse CDF: F^{-1}(u) = (u * range_a + lo_a)^{1/a} -------
        x_right = np.power(u_right * range_a + lo_a, 1.0 / a)
        x_left  = np.power(u_left  * range_a + lo_a, 1.0 / a)

        return np.maximum(np.abs(x_right - values), np.abs(x_left - values))
