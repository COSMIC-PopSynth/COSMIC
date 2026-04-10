"""Vectorized Gaussian Mixture Model for adaptive importance sampling.

Stores all mixture components as numpy arrays and provides vectorized
sampling, PDF evaluation, and EM updates.
"""
import numpy as np
from scipy.stats import multivariate_normal, entropy as scipy_entropy
from .constants import MIN_ENTROPY_CHANGE


class GaussianMixture:
    """A mixture of K multivariate Gaussians.

    Parameters
    ----------
    means : `numpy.ndarray`
        (K, D) array of component means.
    covariances : `numpy.ndarray`
        (K, D, D) array of covariance matrices.
    alphas : `numpy.ndarray`
        (K,) array of mixture weights (must sum to 1).
    rejection_rate : `float`, optional
        Fraction of samples that fall outside bounds, by default 0.0
    """

    def __init__(self, means, covariances, alphas, rejection_rate=0.0):
        self.means = np.asarray(means)
        self.covariances = np.asarray(covariances)
        self.alphas = np.asarray(alphas, dtype=float)
        self.rejection_rate = rejection_rate

    @property
    def n_components(self):
        return len(self.means)

    @property
    def ndim(self):
        return self.means.shape[1]

    @classmethod
    def from_hits(cls, hit_samples, param_space, average_density_one_dim, kappa=1.0):
        """Create a Gaussian mixture by placing one component at each hit.

        Parameters
        ----------
        hit_samples : `numpy.ndarray`
            (K, D) array of hit locations in sampling space.
        param_space : `ParameterSpace`
            Parameter space instance providing bounds and sigma computation.
        average_density_one_dim : `float`
            Characteristic inter-sample spacing,
            ``1 / num_explored ** (1 / D)``.
        kappa : `float`, optional
            Width scaling factor for the Gaussian covariances, by default 1.0

        Returns
        -------
        `GaussianMixture`
            A new mixture with one component centred on each hit.
        """
        K = len(hit_samples)
        D = param_space.ndim

        # Compute sigma for each hit and dimension
        sigmas = param_space.compute_sigma(hit_samples, average_density_one_dim)

        # Build diagonal covariance matrices, scaled by kappa
        covariances = np.zeros((K, D, D))
        for k in range(K):
            covariances[k] = np.diag(sigmas[k] ** 2) * kappa ** 2

        # Equal mixture weights
        alphas = np.full(K, 1.0 / K)

        return cls(hit_samples.copy(), covariances, alphas)

    def sample(self, n_total, param_space, consider_rejection=False, rng=None):
        """Sample from the mixture distribution.

        Parameters
        ----------
        n_total : `int`
            Total number of samples desired.
        param_space : `ParameterSpace`
            Parameter space used for bounds checking.
        consider_rejection : `bool`, optional
            If True, oversample to account for the current rejection
            rate, by default False
        rng : `numpy.random.Generator`, optional
            Random number generator, by default None

        Returns
        -------
        samples : `numpy.ndarray`
            (M, D) array of samples in sampling space.
        mask : `numpy.ndarray`
            (M,) boolean array indicating which samples are in bounds.
        gaussian_indices : `numpy.ndarray`
            (M,) integer array indicating which component generated each
            sample.
        """
        rng = rng or np.random.default_rng()
        all_samples = []
        all_indices = []

        for k in range(self.n_components):
            n_k = int(np.ceil(n_total * self.alphas[k]))
            if consider_rejection and self.rejection_rate < 1.0:
                n_k = int(2 * np.ceil(n_k / (1 - self.rejection_rate)))
            if n_k <= 0:
                continue
            s = rng.multivariate_normal(self.means[k], self.covariances[k], size=n_k)
            all_samples.append(s)
            all_indices.append(np.full(n_k, k, dtype=int))

        if not all_samples:
            D = param_space.ndim
            return np.empty((0, D)), np.empty(0, dtype=bool), np.empty(0, dtype=int)

        samples = np.vstack(all_samples)
        gaussian_indices = np.concatenate(all_indices)
        mask = param_space.in_bounds(samples)

        return samples, mask, gaussian_indices

    def pdf(self, samples):
        """Evaluate the mixture PDF at the given samples.

        Parameters
        ----------
        samples : `numpy.ndarray`
            (N, D) array in sampling space.

        Returns
        -------
        `numpy.ndarray`
            (N,) array of PDF values (weighted sum of components).
        """
        N = len(samples)
        result = np.zeros(N)
        for k in range(self.n_components):
            result += self.alphas[k] * multivariate_normal.pdf(
                samples, self.means[k], self.covariances[k], allow_singular=True
            )
        return result

    def component_pdfs(self, samples):
        """Evaluate each component's PDF at the given samples.

        Parameters
        ----------
        samples : `numpy.ndarray`
            (N, D) array in sampling space.

        Returns
        -------
        `numpy.ndarray`
            (K, N) array where element [k, n] is the PDF of component k
            evaluated at sample n.
        """
        K = self.n_components
        N = len(samples)
        xPDF = np.empty((K, N))
        for k in range(K):
            xPDF[k, :] = multivariate_normal.pdf(
                samples, self.means[k], self.covariances[k], allow_singular=True
            )
        return xPDF

    def compute_rejection_rate(self, param_space, compute_derived_fn, reject_fn,
                               n_per_component=10000, rng=None):
        """Estimate the rejection rate of the mixture.

        Samples from each component, transforms to physical space, applies
        rejection criteria, and computes the weighted rejection rate.

        Parameters
        ----------
        param_space : `ParameterSpace`
            Parameter space for bounds checking and coordinate transforms.
        compute_derived_fn : `callable`
            Function with signature
            ``(samples_physical, param_names) -> dict`` that computes
            derived quantities.
        reject_fn : `callable`
            Function with signature
            ``(samples_physical, derived, param_names) -> bool_mask``
            returning True for rejected systems.
        n_per_component : `int`, optional
            Number of samples per component for the estimate, by default
            10000
        rng : `numpy.random.Generator`, optional
            Random number generator, by default None

        Returns
        -------
        `float`
            Estimated rejection rate (also stored as
            ``self.rejection_rate``).
        """
        rng = rng or np.random.default_rng()
        fractional_rejected = 0.0

        for k in range(self.n_components):
            n = n_per_component
            s = rng.multivariate_normal(self.means[k], self.covariances[k], size=n)

            # Bounds rejection
            bounds_mask = param_space.in_bounds(s)
            rejected = n - np.sum(bounds_mask)

            # Physical rejection on in-bounds samples
            s_valid = s[bounds_mask]
            if len(s_valid) > 0:
                s_physical = param_space.to_physical(s_valid)
                derived = compute_derived_fn(s_physical, param_space.names)
                phys_rejected = reject_fn(s_physical, derived, param_space.names)
                rejected += np.sum(phys_rejected)

            fractional_rejected += rejected * self.alphas[k] / n

        self.rejection_rate = fractional_rejected
        return self.rejection_rate

    def update_em(self, samples, is_hit, prior_probs, prior_fraction_rejected,
                  tolerance=1e-10, entropies=None):
        """Perform one EM-like update of the mixture parameters.

        Parameters
        ----------
        samples : `numpy.ndarray`
            (N, D) array of samples in sampling space.
        is_hit : `numpy.ndarray`
            (N,) boolean or integer array indicating hits.
        prior_probs : `numpy.ndarray`
            (N,) array of prior probabilities for each sample.
        prior_fraction_rejected : `float`
            Fraction of prior samples that are physically rejected.
        tolerance : `float`, optional
            Minimum mixture weight to keep a component, by default 1e-10
        entropies : `list`, optional
            List of previous entropy values (mutated in place for
            convergence tracking), by default None

        Returns
        -------
        `bool`
            True if the entropy check triggers reversion to the previous
            mixture state.
        """
        pi_norm = 1.0 / (1 - prior_fraction_rejected)
        q_norm = 1.0 / (1 - self.rejection_rate)
        pi = prior_probs * pi_norm
        is_hit = np.asarray(is_hit, dtype=float)

        N = len(samples)
        K = self.n_components

        # Compute component PDFs: (K, N)
        xPDF = self.component_pdfs(samples).T  # -> (N, K)
        qPDF = xPDF * self.alphas * q_norm  # (N, K)

        # Responsibilities
        qPDF_sum = np.sum(qPDF, axis=1)  # (N,)
        rho = qPDF / qPDF_sum[:, None]  # (N, K)

        # Importance weights
        gaussian_weights = (pi * is_hit) / qPDF_sum  # (N,)
        weight_sum = np.sum(gaussian_weights)
        if weight_sum == 0:
            return False
        weights_normalized = (gaussian_weights / weight_sum)[:, None]  # (N, 1)

        # Update alphas
        new_alphas = np.sum(weights_normalized * rho, axis=0)  # (K,)

        # Remove insignificant components
        keep = new_alphas > tolerance
        if not np.any(keep):
            return False

        new_alphas = new_alphas[keep]
        rho = rho[:, keep]
        K_new = int(np.sum(keep))

        # Update means
        new_means = np.empty((K_new, self.ndim))
        for d in range(self.ndim):
            new_means[:, d] = np.sum(
                weights_normalized * samples[:, d:d+1] * rho, axis=0
            )
        new_means = new_means / new_alphas[:, None]

        # Update covariances
        old_covs = self.covariances[keep]
        new_covs = np.empty_like(old_covs)
        for k in range(K_new):
            distance = (new_means[k] - samples)[:, :, None]  # (N, D, 1)
            matrix = np.einsum('nij,nji->nij', distance, distance)  # (N, D, D)
            factor = weights_normalized[:, 0] * rho[:, k]  # (N,)
            new_covs[k] = np.sum(factor[:, None, None] * matrix, axis=0) / new_alphas[k]

        # Entropy check
        entropy_change = np.exp(scipy_entropy(weights_normalized[:, 0])) / N
        if entropies is not None:
            if len(entropies) >= 1 and entropy_change - entropies[-1] < MIN_ENTROPY_CHANGE:
                return True  # Signal to revert
            entropies.append(entropy_change)

        # Apply updates
        self.means = new_means
        self.covariances = new_covs
        self.alphas = new_alphas

        return False
